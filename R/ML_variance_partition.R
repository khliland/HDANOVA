.ML_variance_partition_list <- function(models, type = c(1, 2, 3)) {
  stopifnot(all(sapply(models, inherits, what = "merMod")))

  type <- match.arg(type)

  # Force REML = FALSE
  models <- lapply(models, function(model) {
    if (lme4::getME(model, "devcomp")$dims["REML"] != 0){
      suppressMessages(suppressWarnings(lmer(formula(model), data=model@frame, REML=FALSE)))
#      suppressMessages(suppressWarnings(update(model, REML = FALSE)))
    } else model
  })

  # Refit helper
  safe_refit <- function(rhs) {
      suppressMessages(suppressWarnings(lmer(formula(model), data=model@frame, REML=FALSE)))
    fml <- reformulate(rhs, response = response_var)
#    fml <- as.formula(paste(response_var, "~", rhs))
    suppressMessages(suppressWarnings(lmer(fml, data = data, REML = FALSE)))
  }

  for(model in models){
    # Extract design info
    data <- model@frame
    response_var <- all.vars(formula(model))[1]
    full_formula <- formula(model)
    fixed_formula <- lme4::nobars(full_formula)
    fixed_terms <- attr(terms(fixed_formula), "term.labels")

    # Random parts
    random_terms <- lme4::findbars(full_formula)
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
      theta <- fixef(model_full)

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
    vc <- as.data.frame(VarCorr(model_full))
    random_effects <- vc[!is.na(vc$var1), ]
    random_labels <- paste(random_effects$grp, random_effects$var1, sep = ":")
    random_vals <- setNames(random_effects$vcov, random_labels)

    # Residual
    resid_var <- attr(VarCorr(model_full), "sc")^2

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
    suppressMessages(suppressWarnings(lmer(formula(model), data=model@frame, REML=FALSE)))
    # suppressMessages(suppressWarnings(update(model, REML = FALSE)))
  } else model

  # Refit helper
  safe_refit <- function(rhs) {
    fml <- reformulate(rhs, response = response_var)
#    fml <- as.formula(paste(response_var, "~", rhs))
    suppressMessages(suppressWarnings(lmer(fml, data = data, REML = FALSE)))
  }

  # Extract design info
  data <- model@frame
  response_var <- all.vars(formula(model))[1]
  full_formula <- formula(model)
  fixed_formula <- lme4::nobars(full_formula)
  fixed_terms <- attr(terms(fixed_formula), "term.labels")

  # Random parts
  random_terms <- lme4::findbars(full_formula)
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
    theta <- fixef(model_full)

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
  vc <- as.data.frame(VarCorr(model_full))
  random_effects <- vc[!is.na(vc$var1), ]
  random_labels <- paste(random_effects$grp, random_effects$var1, sep = ":")
  random_vals <- setNames(random_effects$vcov, random_labels)

  # Residual
  resid_var <- attr(VarCorr(model_full), "sc")^2

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
