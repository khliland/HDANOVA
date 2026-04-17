.ML_variance_partition_list <- function(models, type = c(1, 2, 3)) {
  stopifnot(all(sapply(models, inherits, what = "merMod")))

  type <- match.arg(type)

  # Force REML = FALSE
  models <- lapply(models, function(model) {
    if (lme4::getME(model, "devcomp")$dims["REML"] != 0){
      suppressMessages(suppressWarnings(lme4::lmer(formula(model), data=model@frame, REML=FALSE)))
#      suppressMessages(suppressWarnings(update(model, REML = FALSE)))
    } else model
  })

  # Refit helper
  safe_refit <- function(rhs) {
      suppressMessages(suppressWarnings(lme4::lmer(formula(model), data=model@frame, REML=FALSE)))
    fml <- reformulate(rhs, response = response_var)
#    fml <- as.formula(paste(response_var, "~", rhs))
    suppressMessages(suppressWarnings(lme4::lmer(fml, data = data, REML = FALSE)))
  }

  for(model in models){
    # Extract design info
    data <- model@frame
    response_var <- all.vars(formula(model))[1]
    full_formula <- formula(model)
    fixed_formula <- reformulas::nobars(full_formula)
    fixed_terms <- attr(terms(fixed_formula), "term.labels")

    # Random parts
    random_terms <- reformulas::findbars(full_formula)
    random_strs <- vapply(random_terms, function(x) paste0("(", deparse(x), ")"), character(1))
    random_rhs <- paste(random_strs, collapse = " + ")

    # Refit full model to be safe
    full_rhs <- paste(c(fixed_terms, random_rhs), collapse = " + ")
    model_full <- safe_refit(full_rhs)
    n <- nrow(data)
    ll_full <- logLik(model_full)
    fitted_vals <- fitted(model_full)
    fixed_var_total <- var(fitted_vals)

    # Null model (random only)
    null_rhs <- paste("1", random_rhs, sep = " + ")
    model_null <- safe_refit(null_rhs)
    ll_null <- logLik(model_null)

    ## --- FIXED EFFECTS ---
    fixed_contrib <- numeric(length(fixed_terms))
    names(fixed_contrib) <- fixed_terms

    if (type %in% c(1, 2)) {
      if (type == 2) {
        for (term in fixed_terms) {
          reduced <- setdiff(fixed_terms, term)
          rhs <- if (length(reduced)) paste(c(reduced, random_rhs), collapse = " + ") else random_rhs
          reduced_model <- try(safe_refit(rhs), silent = TRUE)
          if (inherits(reduced_model, "try-error")) {
            warning(sprintf("Failed to fit model without '%s'", term))
            fixed_contrib[term] <- NA
            next
          }
          ll_reduced <- logLik(reduced_model)
          ll_gain <- as.numeric(ll_full - ll_reduced)
          fixed_contrib[term] <- (ll_gain / as.numeric(ll_full - ll_null)) * fixed_var_total
        }
      }

      if (type == 1) {
        current_terms <- character(0)
        ll_prev <- ll_null
        for (term in fixed_terms) {
          current_terms <- c(current_terms, term)
          rhs <- paste(c(current_terms, random_rhs), collapse = " + ")
          step_model <- try(safe_refit(rhs), silent = TRUE)
          if (inherits(step_model, "try-error")) {
            warning(sprintf("Failed to fit model adding '%s'", term))
            fixed_contrib[term] <- NA
            next
          }
          ll_step <- logLik(step_model)
          ll_gain <- as.numeric(ll_step - ll_prev)
          ll_prev <- ll_step
          fixed_contrib[term] <- (ll_gain / as.numeric(ll_full - ll_null)) * fixed_var_total
        }
      }

    } else if (type == 3) {
      # Use L-matrix method
      X <- model.matrix(model_full)
      assign_vec <- attr(X, "assign")
      term_labels <- attr(terms(model_full), "term.labels")
      has_intercept <- any(assign_vec == 0)
      if (has_intercept) term_labels <- c("(Intercept)", term_labels)

      XtX_inv <- MASS::ginv(t(X) %*% X)
      theta <- lme4::fixef(model_full)

      get_L <- function(term_idx) {
        cols <- which(assign_vec == term_idx)
        if (length(cols) == 0) return(matrix(0, 0, ncol(X)))
        L <- matrix(0, nrow = length(cols), ncol = ncol(X))
        for (i in seq_along(cols)) L[i, cols[i]] <- 1
        rownames(L) <- colnames(X)[cols]
        colnames(L) <- colnames(X)
        return(L)
      }

      ss_list <- numeric(length(term_labels))
      names(ss_list) <- term_labels
      for (i in seq_along(term_labels)) {
        term_idx <- ifelse(has_intercept, i - 1, i)
        L <- get_L(term_idx)
        if (nrow(L) == 0) {
          ss_list[i] <- 0
          next
        }
        M <- L %*% XtX_inv %*% t(L)
        M_inv <- MASS::ginv(M)
        ss_list[i] <- as.numeric((t(theta) %*% t(L) %*% M_inv %*% L %*% theta) / n)
      }

      fixed_contrib <- ss_list
      if (has_intercept) {
        fixed_contrib <- fixed_contrib[names(fixed_contrib) != "(Intercept)"]
      }
    }

    ## --- RANDOM EFFECTS ---
    vc <- as.data.frame(lme4::VarCorr(model_full))
    random_effects <- vc[!is.na(vc$var1), ]
    random_labels <- paste(random_effects$grp, random_effects$var1, sep = ":")
    random_vals <- setNames(random_effects$vcov, random_labels)

    # Residual
    resid_var <- attr(lme4::VarCorr(model_full), "sc")^2

    ## --- COMBINE ---
    if(!exists("all_parts")){
      all_parts <- c(fixed_contrib, random_vals, Residual = resid_var)
    } else {
      all_parts <- all_parts + c(fixed_contrib, random_vals, Residual = resid_var)
    }
  }
  ## Wrap-up
  total <- sum(all_parts, na.rm = TRUE)
  proportions <- all_parts / total

  result <- data.frame(
    Term = names(all_parts),
    Type = c(rep("Fixed", length(fixed_contrib)),
             rep("Random", length(random_vals)),
             "Residual"),
    Variance = as.numeric(all_parts),
    Proportion = proportions
  )

  rownames(result) <- NULL
  return(result)
}

.ML_variance_partition_single <- function(model, type = c(1, 2, 3)) {
  stopifnot(inherits(model, "merMod"))

  # Force REML = FALSE
  if (lme4::getME(model, "devcomp")$dims["REML"] != 0){
    suppressMessages(suppressWarnings(lme4::lmer(formula(model), data=model@frame, REML=FALSE)))
    # suppressMessages(suppressWarnings(update(model, REML = FALSE)))
  } else model

  # Refit helper
  safe_refit <- function(rhs) {
    fml <- reformulate(rhs, response = response_var)
#    fml <- as.formula(paste(response_var, "~", rhs))
    suppressMessages(suppressWarnings(lme4::lmer(fml, data = data, REML = FALSE)))
  }

  # Extract design info
  data <- model@frame
  response_var <- all.vars(formula(model))[1]
  full_formula <- formula(model)
  fixed_formula <- reformulas::nobars(full_formula)
  fixed_terms <- attr(terms(fixed_formula), "term.labels")

  # Random parts
  random_terms <- reformulas::findbars(full_formula)
  random_strs <- vapply(random_terms, function(x) paste0("(", deparse(x), ")"), character(1))
  random_rhs <- paste(random_strs, collapse = " + ")

  # Refit full model to be safe
  full_rhs <- paste(c(fixed_terms, random_rhs), collapse = " + ")
  model_full <- safe_refit(full_rhs)
  n <- nrow(data)
  ll_full <- logLik(model_full)
  fitted_vals <- fitted(model_full)
  fixed_var_total <- var(fitted_vals)

  # Null model (random only)
  null_rhs <- paste("1", random_rhs, sep = " + ")
  model_null <- safe_refit(null_rhs)
  ll_null <- logLik(model_null)

  ## --- FIXED EFFECTS ---
  fixed_contrib <- numeric(length(fixed_terms))
  names(fixed_contrib) <- fixed_terms

  if (type %in% c(1, 2)) {
    if (type == 2) {
      for (term in fixed_terms) {
        reduced <- setdiff(fixed_terms, term)
        rhs <- if (length(reduced)) paste(c(reduced, random_rhs), collapse = " + ") else random_rhs
        reduced_model <- try(safe_refit(rhs), silent = TRUE)
        if (inherits(reduced_model, "try-error")) {
          warning(sprintf("Failed to fit model without '%s'", term))
          fixed_contrib[term] <- NA
          next
        }
        ll_reduced <- logLik(reduced_model)
        ll_gain <- as.numeric(ll_full - ll_reduced)
        fixed_contrib[term] <- (ll_gain / as.numeric(ll_full - ll_null)) * fixed_var_total
      }
    }

    if (type == 1) {
      current_terms <- character(0)
      ll_prev <- ll_null
      for (term in fixed_terms) {
        current_terms <- c(current_terms, term)
        rhs <- paste(c(current_terms, random_rhs), collapse = " + ")
        step_model <- try(safe_refit(rhs), silent = TRUE)
        if (inherits(step_model, "try-error")) {
          warning(sprintf("Failed to fit model adding '%s'", term))
          fixed_contrib[term] <- NA
          next
        }
        ll_step <- logLik(step_model)
        ll_gain <- as.numeric(ll_step - ll_prev)
        ll_prev <- ll_step
        fixed_contrib[term] <- (ll_gain / as.numeric(ll_full - ll_null)) * fixed_var_total
      }
    }

  } else if (type == 3) {
    # Use L-matrix method
    X <- model.matrix(model_full)
    assign_vec <- attr(X, "assign")
    term_labels <- attr(terms(model_full), "term.labels")
    has_intercept <- any(assign_vec == 0)
    if (has_intercept) term_labels <- c("(Intercept)", term_labels)

    XtX_inv <- MASS::ginv(t(X) %*% X)
    theta <- lme4::fixef(model_full)

    get_L <- function(term_idx) {
      cols <- which(assign_vec == term_idx)
      if (length(cols) == 0) return(matrix(0, 0, ncol(X)))
      L <- matrix(0, nrow = length(cols), ncol = ncol(X))
      for (i in seq_along(cols)) L[i, cols[i]] <- 1
      rownames(L) <- colnames(X)[cols]
      colnames(L) <- colnames(X)
      return(L)
    }

    ss_list <- numeric(length(term_labels))
    names(ss_list) <- term_labels
    for (i in seq_along(term_labels)) {
      term_idx <- ifelse(has_intercept, i - 1, i)
      L <- get_L(term_idx)
      if (nrow(L) == 0) {
        ss_list[i] <- 0
        next
      }
      M <- L %*% XtX_inv %*% t(L)
      M_inv <- MASS::ginv(M)
      ss_list[i] <- as.numeric((t(theta) %*% t(L) %*% M_inv %*% L %*% theta) / n)
    }

    fixed_contrib <- ss_list
    if (has_intercept) {
      fixed_contrib <- fixed_contrib[names(fixed_contrib) != "(Intercept)"]
    }
  }

  ## --- RANDOM EFFECTS ---
  vc <- as.data.frame(lme4::VarCorr(model_full))
  random_effects <- vc[!is.na(vc$var1), ]
  random_labels <- paste(random_effects$grp, random_effects$var1, sep = ":")
  random_vals <- setNames(random_effects$vcov, random_labels)

  # Residual
  resid_var <- attr(lme4::VarCorr(model_full), "sc")^2

  ## --- COMBINE ---
  all_parts <- c(fixed_contrib, random_vals, Residual = resid_var)
  total <- sum(all_parts, na.rm = TRUE)
  proportions <- all_parts / total

  result <- data.frame(
    Term = names(all_parts),
    Type = c(rep("Fixed", length(fixed_contrib)),
             rep("Random", length(random_vals)),
             "Residual"),
    Variance = as.numeric(all_parts),
    Proportion = proportions
  )

  rownames(result) <- NULL
  return(result)
}

.ML_variance_partition_all_types <- function(models,
                                             type = 3,
                                             method = c("exact_refit", "wald", "ls"),
                                             control = lme4::lmerControl(),
                                             verbose = FALSE) {
  method_label <- match.arg(method, c("exact_refit", "wald", "ls"))
  method_kernel <- method_label
  type <- .hda_normalize_ss_type(type)
  models <- .hda_normalize_models_input(models)
  .hda_validate_model_family(models)

  template_info <- .hda_prepare_refit_templates(models[[1]], type)
  all_results <- vector("list", length(models))

  for(m_idx in seq_along(models)){
    m <- models[[m_idx]]
    n_obs <- stats::nobs(m)
    fixed_terms <- template_info$fixed_terms
    fixed_ss <- setNames(numeric(length(fixed_terms)), fixed_terms)

    if(method_kernel == "exact_refit"){
      fixed_ss <- .hda_fixed_ss_exact_refit(
        model = m,
        template_info = template_info,
        control = control,
        verbose = verbose
      )
    } else if(method_kernel == "wald"){
      fixed_ss <- .hda_fixed_ss_wald(
        model = m,
        type = type,
        fixed_terms = fixed_terms
      )
    } else {
      fixed_ss <- .hda_fixed_ss_ls(model = m, fixed_terms = fixed_terms)
    }

    random_resid <- .hda_random_residual_ss(model = m, n_obs = n_obs)
    all_results[[m_idx]] <- c(fixed_ss, random_resid)
  }

  all_names <- unique(unlist(lapply(all_results, names)))
  ss_matrix <- matrix(0, nrow = length(models), ncol = length(all_names),
                      dimnames = list(NULL, all_names))
  for(i in seq_along(all_results)) {
    match_idx <- match(names(all_results[[i]]), all_names)
    ss_matrix[i, match_idx] <- all_results[[i]]
  }

  total_ssq <- colSums(ss_matrix, na.rm = TRUE)
  total_ss <- sum(total_ssq, na.rm = TRUE)
  term_type <- .hda_partition_term_type(all_names, template_info$fixed_terms)
  proportions <- if(total_ss > 0) as.numeric(total_ssq / total_ss) else rep(NA_real_, length(total_ssq))

  data.frame(
    Term = names(total_ssq),
    Type = term_type,
    SSQ = as.numeric(total_ssq),
    Proportion = proportions,
    Method = rep(method_label, length(total_ssq)),
    SStype = rep(type, length(total_ssq)),
    stringsAsFactors = FALSE
  )
}


.hda_normalize_ss_type <- function(type){
  if(length(type) != 1 || is.na(type))
    stop("'type' must be a single value in {1, 2, 3}.")
  if(is.character(type))
    type <- as.integer(type)
  if(!is.numeric(type) || !(type %in% c(1, 2, 3)))
    stop("'type' must be 1, 2, or 3.")
  as.integer(type)
}


.hda_normalize_models_input <- function(models){
  if(inherits(models, "merMod"))
    models <- list(models)
  if(!is.list(models) || length(models) == 0)
    stop("'models' must be a non-empty list of merMod objects (or a single merMod).")
  if(!all(vapply(models, inherits, what = "merMod", FUN.VALUE = logical(1))))
    stop("All entries in 'models' must inherit from 'merMod'.")
  models
}


.hda_validate_model_family <- function(models){
  fixed_ref <- attr(terms(models[[1]]), "term.labels")
  rand_ref <- vapply(reformulas::findbars(formula(models[[1]])), deparse, character(1))
  for(i in seq_along(models)){
    fixed_i <- attr(terms(models[[i]]), "term.labels")
    rand_i <- vapply(reformulas::findbars(formula(models[[i]])), deparse, character(1))
    if(!identical(fixed_i, fixed_ref) || !identical(rand_i, rand_ref)){
      stop("All models must share the same fixed and random structure for aggregation.")
    }
  }
  invisible(TRUE)
}


.hda_prepare_refit_templates <- function(model, type){
  full_formula <- formula(model)
  response_var <- all.vars(full_formula)[1]
  fixed_formula <- reformulas::nobars(full_formula)
  fixed_terms <- attr(terms(fixed_formula), "term.labels")
  random_terms <- reformulas::findbars(full_formula)
  random_strs <- vapply(random_terms, function(x) paste0("(", deparse(x), ")"), character(1))
  random_rhs <- if(length(random_strs) > 0) paste(random_strs, collapse = " + ") else ""
  data_ref <- model@frame

  build_formula <- function(terms_fixed){
    rhs_parts <- c(terms_fixed, random_rhs)
    rhs_parts <- rhs_parts[rhs_parts != ""]
    rhs <- if(length(rhs_parts) > 0) paste(rhs_parts, collapse = " + ") else "1"
    stats::reformulate(rhs, response = response_var)
  }

  # Build per-term full/reduced fixed-term sets according to SS type.
  factors <- attr(terms(fixed_formula), "factors")
  term_sets <- vector("list", length(fixed_terms))
  names(term_sets) <- fixed_terms
  for(idx in seq_along(fixed_terms)){
    if(type == 1){
      full_terms <- fixed_terms[seq_len(idx)]
      reduced_terms <- if(idx == 1) character(0) else fixed_terms[seq_len(idx - 1)]
    } else if(type == 2){
      current_comp <- factors[, idx]
      is_descendant <- apply(factors, 2, function(col) all(current_comp <= col) && any(current_comp < col))
      descendants <- which(is_descendant)
      full_terms <- fixed_terms[setdiff(seq_along(fixed_terms), descendants)]
      reduced_terms <- setdiff(full_terms, fixed_terms[idx])
    } else {
      full_terms <- fixed_terms
      reduced_terms <- setdiff(fixed_terms, fixed_terms[idx])
    }
    term_sets[[idx]] <- list(full_terms = full_terms, reduced_terms = reduced_terms)
  }

  full_formula_global <- build_formula(fixed_terms)
  null_formula <- build_formula(character(0))
  template_full <- suppressMessages(suppressWarnings(lme4::lmer(full_formula_global, data = data_ref, REML = FALSE)))
  template_null <- suppressMessages(suppressWarnings(lme4::lmer(null_formula, data = data_ref, REML = FALSE)))

  term_templates <- vector("list", length(fixed_terms))
  names(term_templates) <- fixed_terms
  for(idx in seq_along(fixed_terms)){
    fs <- term_sets[[idx]]
    form_full <- build_formula(fs$full_terms)
    form_red <- build_formula(fs$reduced_terms)
    term_templates[[idx]] <- list(
      full = suppressMessages(suppressWarnings(lme4::lmer(form_full, data = data_ref, REML = FALSE))),
      reduced = suppressMessages(suppressWarnings(lme4::lmer(form_red, data = data_ref, REML = FALSE)))
    )
  }

  # Pre-extract design matrix info once (X and assign_vec are response-invariant)
  X_ref <- lme4::getME(template_full, "X")
  assign_vec_ref <- attr(X_ref, "assign")

  list(
    fixed_terms = fixed_terms,
    response_var = response_var,
    template_full = template_full,
    template_null = template_null,
    term_templates = term_templates,
    X_ref = X_ref,
    assign_vec_ref = assign_vec_ref,
    type = type
  )
}


.hda_refit_or_update <- function(template, y, control, verbose){
  refit_res <- try(suppressMessages(suppressWarnings(lme4::refit(template, newresp = y))), silent = TRUE)
  if(!inherits(refit_res, "try-error"))
    return(refit_res)
  if(isTRUE(verbose))
    warning("refit() failed; falling back to fresh lmer() fit for one response.")
  form <- formula(template)
  data_ref <- template@frame
  data_ref[[all.vars(form)[1]]] <- y
  suppressMessages(suppressWarnings(lme4::lmer(form, data = data_ref, REML = FALSE, control = control)))
}


.hda_fixed_ss_exact_refit <- function(model, template_info, control, verbose){
  y <- model@frame[[all.vars(formula(model))[1]]]
  fixed_terms <- template_info$fixed_terms
  contrib <- setNames(numeric(length(fixed_terms)), fixed_terms)
  type <- template_info$type
  X <- template_info$X_ref           # design matrix is response-invariant
  assign_vec <- template_info$assign_vec_ref

  if(type == 3L){
    # Type III: LS-norm approach — avoids non-hierarchical reduced-model artefacts.
    # SS_k = ||X_k * beta_k||^2  (marginal contribution of term k in the full model)
    full_m <- .hda_refit_or_update(template_info$template_full, y, control, verbose)
    beta <- lme4::fixef(full_m)
    for(idx in seq_along(fixed_terms)){
      cols <- which(assign_vec == idx)
      if(length(cols) == 0){ contrib[idx] <- 0; next }
      contrib[idx] <- sum((X[, cols, drop = FALSE] %*% beta[cols])^2)
    }
  } else {
    # Type I / II: hierarchically nested model comparison — fitted-value differencing.
    for(idx in seq_along(fixed_terms)){
      tm <- template_info$term_templates[[idx]]
      fit_full <- .hda_refit_or_update(tm$full, y, control, verbose)
      fit_red  <- .hda_refit_or_update(tm$reduced, y, control, verbose)
      contrib[idx] <- sum((stats::fitted(fit_full) - stats::fitted(fit_red))^2)
    }
  }
  contrib
}


.hda_safe_inverse <- function(mat){
  inv <- try(solve(mat), silent = TRUE)
  if(inherits(inv, "try-error"))
    inv <- MASS::ginv(mat)
  inv
}


.hda_fixed_ss_wald <- function(model, type, fixed_terms){
  X <- lme4::getME(model, "X")
  beta_full <- lme4::fixef(model)
  assign_vec <- attr(X, "assign")
  RX <- lme4::getME(model, "RX")
  inf_full <- t(RX) %*% RX
  contrib <- setNames(numeric(length(fixed_terms)), fixed_terms)

  if(type == 1){
    L_chol <- t(chol(inf_full))
    seq_components <- (t(L_chol) %*% beta_full)^2
    for(idx in seq_along(fixed_terms))
      contrib[idx] <- sum(seq_components[which(assign_vec == idx)])
    return(contrib)
  }

  factors <- attr(terms(model), "factors")
  cov_full <- .hda_safe_inverse(inf_full)
  for(idx in seq_along(fixed_terms)){
    if(type == 2){
      current_comp <- factors[, idx]
      is_descendant <- apply(factors, 2, function(col) all(current_comp <= col) && any(current_comp < col))
      descendants <- which(is_descendant)
      keep_cols <- which(!(assign_vec %in% descendants))
      cov_sub <- .hda_safe_inverse(inf_full[keep_cols, keep_cols, drop = FALSE])
      beta_sub <- beta_full[keep_cols]
      assign_sub <- assign_vec[keep_cols]
      L_cols <- which(assign_sub == idx)
      L <- matrix(0, nrow = length(L_cols), ncol = length(keep_cols))
      for(j in seq_along(L_cols))
        L[j, L_cols[j]] <- 1
      middle <- L %*% cov_sub %*% t(L)
      contrib[idx] <- as.numeric(t(L %*% beta_sub) %*% .hda_safe_inverse(middle) %*% (L %*% beta_sub))
    } else {
      cols <- which(assign_vec == idx)
      L <- matrix(0, nrow = length(cols), ncol = ncol(X))
      for(j in seq_along(cols))
        L[j, cols[j]] <- 1
      middle <- L %*% cov_full %*% t(L)
      contrib[idx] <- as.numeric(t(L %*% beta_full) %*% .hda_safe_inverse(middle) %*% (L %*% beta_full))
    }
  }
  contrib
}


.hda_fixed_ss_ls <- function(model, fixed_terms){
  X <- lme4::getME(model, "X")
  beta <- lme4::fixef(model)
  assign_vec <- attr(X, "assign")
  contrib <- setNames(numeric(length(fixed_terms)), fixed_terms)
  for(idx in seq_along(fixed_terms)){
    cols <- which(assign_vec == idx)
    if(length(cols) == 0){
      contrib[idx] <- 0
    } else {
      fitted_term <- X[, cols, drop = FALSE] %*% beta[cols]
      contrib[idx] <- sum(fitted_term^2)
    }
  }
  contrib
}


.hda_random_residual_ss <- function(model, n_obs){
  # Exact residual SS
  resid_ss <- c(Residuals = sum(stats::residuals(model)^2))

  # Per-grouping-factor BLUP contribution: sum((Z_g * b_g)^2)
  random_bars <- reformulas::findbars(formula(model))
  if(length(random_bars) == 0) return(resid_ss)

  fixed_pred <- suppressWarnings(stats::predict(model, re.form = NA))

  random_ss <- vapply(random_bars, function(bar){
    re_form <- stats::reformulate(paste0("(", deparse(bar), ")"))
    pred_with <- tryCatch(
      suppressWarnings(stats::predict(model, re.form = re_form)),
      error = function(e) fixed_pred
    )
    sum((pred_with - fixed_pred)^2)
  }, numeric(1))

  # Label = RHS of bar expression (grouping factor name)
  bar_labels <- vapply(random_bars, function(bar){
    parts <- strsplit(deparse(bar), "\\s*\\|\\s*")[[1]]
    trimws(parts[length(parts)])
  }, character(1))

  names(random_ss) <- bar_labels
  c(random_ss, resid_ss)
}


.hda_partition_term_type <- function(term_names, fixed_terms){
  vapply(term_names, function(x){
    if(x %in% c("Residual", "Residuals"))
      return("Residual")
    if(x %in% fixed_terms)
      return("Fixed")
    "Random"
  }, FUN.VALUE = character(1))
}