#' @title High-Dimensional Analysis of Variance (QR Prototype)
#'
#' @description
#' Experimental QR-based alternative to \code{hdanova()}.
#' Phase 1 supports fixed-effect linear models only.
#'
#' @inheritParams hdanova
#' @param respect_SStype Logical; if \code{FALSE} (default), keep the legacy
#' regression-based effect matrices in \code{LS}. If \code{TRUE}, expose the
#' SS-type-aligned QR effect matrices in \code{LS} while retaining the legacy
#' matrices in \code{LS_regression}.
#'
#' @return
#' An \code{hdanova}-compatible object.
#'
#' @export
hdanova2 <- function(formula, data, subset, weights, na.action, family,
                     unrestricted = FALSE,
                     add_error = FALSE,
                     aug_error = "denominator",
                     use_ED = FALSE,
                     pca.in = FALSE,
                     contrasts = "contr.sum",
                     equal_baseline = FALSE,
                     SStype = "II",
                     REML = NULL,
                     scale = FALSE,
                     respect_SStype = FALSE){

  # Simplify SStype in the same style as hdanova()
  if(is.character(SStype))
    SStype <- nchar(SStype)
  if(!SStype %in% 1:3)
    stop("'SStype' must be one of 'I', 'II' or 'III'.")

  # Phase-1 scope guardrails
  if(!missing(family))
    stop("hdanova2() Phase 1 supports LM only. Use hdanova() for GLM/GLMM models.")
  if(any(grepl("|", formula, fixed = TRUE)))
    stop("hdanova2() currently supports mixlm-style random effects r(), not lme4 '|' notation.")
  mixed_r <- any(grepl("r(", formula, fixed = TRUE))
  if(!is.logical(unrestricted) || length(unrestricted) != 1 || is.na(unrestricted))
    stop("'unrestricted' must be TRUE or FALSE.")
  if(!is.logical(respect_SStype) || length(respect_SStype) != 1 || is.na(respect_SStype))
    stop("'respect_SStype' must be TRUE or FALSE.")
  aug_is_residual <- is.character(aug_error) && length(aug_error) == 1 &&
    aug_error %in% c("residual", "residuals")
  aug_is_denominator <- is.character(aug_error) && length(aug_error) == 1 &&
    identical(aug_error, "denominator")
  aug_is_numeric <- is.numeric(aug_error) && length(aug_error) == 1 && is.finite(aug_error)
  if(!(aug_is_residual || aug_is_denominator || aug_is_numeric))
    stop("'aug_error' must be 'denominator', 'residual'/'residuals', or a numeric alpha in [0,1].")
  if(aug_is_numeric && (aug_error < 0 || aug_error > 1))
    stop("Numeric 'aug_error' must be in [0,1].")
  if(!is.logical(use_ED) || length(use_ED) != 1 || is.na(use_ED))
    stop("'use_ED' must be TRUE or FALSE.")
  if(!is.null(REML) && !is.logical(REML))
    stop("'REML' must be NULL, TRUE or FALSE.")
  if(is.logical(REML) && !mixed_r)
    stop("'REML' is only supported for mixed models with r() terms.")
  use_ED_active <- isTRUE(use_ED) && mixed_r && is.logical(REML) && aug_is_numeric
  if(isTRUE(use_ED) && !use_ED_active)
    warning("'use_ED' currently only affects numeric 'aug_error' in REML/ML mixed-model fits; using standard degrees of freedom instead.")

  # Set contrasts in line with hdanova()
  old_options <- options(contrasts = c(contrasts, "contr.poly"))
  on.exit(options(old_options), add = TRUE)

  # Build model.frame once to keep subset/weights/na.action behavior consistent.
  formula_mf <- formula
  if(mixed_r){
    rw <- mixlm::random.worker(formula, data, REML)
    formula_mf <- rw$formula
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$formula <- formula_mf
  mf[[1]] <- as.name("model.frame")
  mf_data <- eval(mf, envir = parent.frame())

  Y <- stats::model.response(mf_data)
  if(inherits(Y, "AsIs"))
    Y <- unclass(Y)
  Yorig <- Y <- as.matrix(Y)
  N <- nrow(Y)
  p <- ncol(Y)
  if(is.null(p))
    stop("Response must be a matrix.")

  ########################### Scale response matrix ##########################
  if(is.logical(scale) && scale){
    # Autoscaling
    Y <- scale(Y, center = TRUE, scale = TRUE)
  } else if(is.numeric(scale)){
    # Scale by numeric value
    if(length(scale) != p)
      stop("Numeric 'scale' must have the same length as the number of variables in Y")
    Y <- scale(Y, center = FALSE, scale = scale)
  } else if(is.character(scale)){
    # Scale by factor and level
    if(length(scale) != 2)
      stop("Character 'scale' must have two elements: factor name and level")
    if(!scale[1] %in% names(mf_data))
      stop(paste0("Factor '", scale[1], "' not found in data"))
    fac <- mf_data[[scale[1]]]
    fac_chr <- as.character(fac)
    if(!scale[2] %in% fac_chr)
      stop(paste0("Level '", scale[2], "' not found in factor '", scale[1], "'"))
    sd <- apply(Y[fac_chr == scale[2],,drop=FALSE], 2, stats::sd)
    Y <- scale(Y, center = FALSE, scale = sd)
  }

  ########################## PCA of response matrix ##########################
  if(pca.in != 0){
    if(is.logical(pca.in) && pca.in)
      pca.in <- which.min(.PCAcv(Y))
    if(pca.in < 1){
      pca <- .pca(Y)
      pca.in <- min(which(cumsum(pca$explvar / 100) >= pca.in))
    }
    if(pca.in > p){
      warning(paste0("Reducing 'pca.in' from ", pca.in, " to the number of variables (", p, ")"))
      pca.in <- p
    }
    pca <- .pca(Y, ncomp = pca.in)
    Y <- pca$scores
  }

  # Replace response with the (possibly scaled) matrix response.
  mf_data[[1]] <- Y
  p_fit <- ncol(Y)

  fit.type <- if(mixed_r) "'lmm' (Linear Mixed Model)" else "'lm' (Linear Model)"

  # Fit model(s) and derive model structures.
  weights_vec <- stats::model.weights(mf_data)
  if(!is.null(weights_vec)){
    weights_vec <- as.numeric(weights_vec)
    if(any(!is.finite(weights_vec)) || any(weights_vec < 0))
      stop("'weights' must be finite and non-negative.")
  }
  models <- NULL
  mom_anova <- NULL
  if(is.logical(REML)){
    resp_name <- as.character(formula[[2]])
    build_reml_args <- function(ymat){
      dat_i <- mf_data
      dat_i[[resp_name]] <- ymat
      args <- list(formula = formula,
                   data = dat_i,
                   REML = REML,
                   contrasts = contrasts)
      if(!is.null(weights_vec))
        args$weights <- weights_vec
      args
    }

    base_mod <- do.call(mixlm::lm, build_reml_args(Y[, 1, drop = FALSE]))
    models <- vector("list", p_fit)
    names(models) <- colnames(Y)

    for(i in seq_len(p_fit)){
      if(i == 1){
        modi <- base_mod
      } else {
        modi <- do.call(mixlm::lm, build_reml_args(Y[, i, drop = FALSE]))
      }
      models[[i]] <- modi
      u <- unlist(lme4::ranef(modi))
      cfi <- c(lme4::fixef(modi), u)
      if(i == 1){
        coefs <- matrix(0.0, length(cfi), p_fit)
        rownames(coefs) <- names(cfi)
      }
      coefs[, i] <- cfi
    }
    mod <- base_mod

    # Denominator mapping for REML/ML follows the legacy path via mixlm least-squares model.
    ls_dat <- mf_data
    ls_dat[[resp_name]] <- Y[, 1, drop = FALSE]
    ls_args <- list(formula = formula,
            data = ls_dat,
            unrestricted = unrestricted,
            equal_baseline = equal_baseline,
            contrasts = contrasts)
    mod_ls <- do.call(mixlm::lm, ls_args)
    mom_anova <- mixlm::AnovaMix(mod_ls, SStype = SStype)
  } else {
    lm_args <- list(formula = formula, data = mf_data, unrestricted = unrestricted,
                    equal_baseline = equal_baseline, contrasts = contrasts)
    if(!is.null(weights_vec))
      lm_args$weights <- weights_vec
    mod <- do.call(mixlm::lm, lm_args)
    if(mixed_r)
      mom_anova <- mixlm::Anova(mod, type = c("I", "II", "III")[SStype])
  }

  M <- model.matrix(mod)
  X_fixed <- M
  assign <- attr(M, "assign")
  modFra <- HDANOVA::extended.model.frame(model.frame(mod), mf_data)
  effs <- attr(terms(mod), "term.labels")

  # Match legacy hdanova behavior: append random-effect design columns for REML/ML branch.
  if(is.logical(REML)){
    Z <- as.matrix(lme4::getME(mod, "Z"))
    M <- cbind(M, Z)
    u <- lme4::ranef(mod)
    effs <- c(effs, names(u))
    max_assign <- max(assign)
    for(i in seq_along(u)){
      max_assign <- max_assign + 1
      assign <- c(assign, rep(max_assign, nrow(u[[i]])))
    }
  }

  n_terms <- length(effs)
  if(n_terms == 0)
    stop("No factors in model")

  # Keep assign/effect mapping untouched for downstream consumers.
  maxDir <- numeric(n_terms)
  for(i in seq_len(n_terms))
    maxDir[i] <- min(sum(assign == i), p)

  qr_res <- .calc_SS_qr_loop(M, Y, assign, effs, SStype, weights = weights_vec)
  SS_matrix <- qr_res$SS_matrix
  Y_hat_full <- qr_res$Y_hat_full
  QR_full <- qr_res$QR_full
  residuals <- Y - Y_hat_full
  obs_weights <- if(is.null(weights_vec)) rep(1, N) else weights_vec
  SSE_full <- colSums(obs_weights * residuals^2)

  df_res <- N - QR_full$rank
  df_terms <- qr_res$df_terms
  MS_matrix <- sweep(SS_matrix, 1, df_terms, "/")
  MS_res <- SSE_full / df_res

  if(!exists("coefs")){
    coefs <- as.matrix(stats::coef(mod))
    if(ncol(coefs) != p_fit)
      coefs <- matrix(coefs, nrow = nrow(coefs), ncol = p_fit)
  }
  colnames(coefs) <- colnames(Y)

  LS_regression <- vector("list", n_terms)
  names(LS_regression) <- effs
  for(i in seq_len(n_terms)){
    estimable <- !is.na(rowSums(coefs[assign == i, , drop = FALSE]))
    LS_regression[[effs[i]]] <- M[, assign == i, drop = FALSE][, estimable, drop = FALSE] %*%
      coefs[assign == i, , drop = FALSE][estimable, , drop = FALSE]
    colnames(LS_regression[[effs[i]]]) <- colnames(Y)
  }
  LS_SStype <- qr_res$effects
  estimable <- !is.na(rowSums(coefs))
  residuals <- Y - M[, estimable, drop = FALSE] %*% coefs[estimable, , drop = FALSE]
  SSE_full <- colSums(obs_weights * residuals^2)
  MS_res <- SSE_full / df_res
  den_structure <- .hda_collect_effect_structure(
    formula = formula,
    term_labels = effs,
    model_frame = modFra,
    mixed_r = mixed_r,
    unrestricted = unrestricted,
    random_info = if(exists("rw")) rw else NULL,
    mom_anova = mom_anova
  )
  den_candidates <- .hda_find_denominator_candidates(den_structure)
  den_rules <- .hda_build_denominator_rules(
    structure = den_structure,
    candidates = den_candidates,
    mom_anova = mom_anova,
    df_res = df_res
  )
  err_map <- .hda_apply_denominator_rules(
    rules = den_rules,
    MS_matrix = MS_matrix,
    MS_res = MS_res,
    df_terms = df_terms,
    df_res = df_res,
    term_labels = effs
  )
  MS_error_matrix <- err_map$MS_error_matrix
  den_dfs <- err_map$den_dfs
  den_labels <- err_map$den_labels

  F_matrix <- P_matrix <- matrix(NA_real_, nrow = n_terms, ncol = p_fit,
                                 dimnames = list(effs, colnames(Y)))
  for(i in seq_len(n_terms)){
    F_matrix[i, ] <- MS_matrix[i, ] / MS_error_matrix[i, ]
    P_matrix[i, ] <- stats::pf(F_matrix[i, ], df_terms[i], den_dfs[i], lower.tail = FALSE)
  }

  # Build per-response ANOVA-like objects to preserve downstream expectations.
  anovas <- vector("list", p_fit)
  names(anovas) <- colnames(Y)
  for(i in seq_len(p_fit)){
    ano <- data.frame("Df" = c(df_terms, df_res),
                      "Sum Sq" = c(SS_matrix[, i], SSE_full[i]),
                      "Mean Sq" = c(MS_matrix[, i], MS_res[i]),
                      "F value" = c(F_matrix[, i], NA_real_),
                      "Pr(>F)" = c(P_matrix[, i], NA_real_),
                      "Error Term" = c(den_labels, NA_character_),
                      check.names = FALSE)
    rownames(ano) <- c(effs, "Residuals")
    anovas[[i]] <- ano
  }

  effects <- setNames(lapply(effs, function(eff) modFra[[eff]]), effs)
  aug_df_terms <- df_terms
  if(mixed_r && !is.null(mom_anova) && !is.null(mom_anova$anova)){
    mom_df <- stats::setNames(as.numeric(mom_anova$anova[["Df"]]), rownames(mom_anova$anova))
    matched <- intersect(names(aug_df_terms), names(mom_df))
    aug_df_terms[matched] <- mom_df[matched]
  }
  EDall <- NULL
  if(use_ED_active){
    ed_df_template <- c(aug_df_terms, Residuals = df_res)
    EDall <- .hda_compute_effective_dimensions(models = models, Y = Y, X = X_fixed, Z = Z,
                                               df_template = ed_df_template)
  }
  aug_regression <- .hda_build_augmented_effect_matrices(
    LS_base = LS_regression,
    residuals = residuals,
    effs = effs,
    mixed_r = mixed_r,
    aug_is_residual = aug_is_residual,
    aug_is_numeric = aug_is_numeric,
    aug_error = aug_error,
    df_terms = aug_df_terms,
    den_dfs = den_dfs,
    den_labels = den_labels,
    den_rules = den_rules,
    use_ED = use_ED_active,
    EDall = EDall
  )
  aug_SStype <- .hda_build_augmented_effect_matrices(
    LS_base = LS_SStype,
    residuals = residuals,
    effs = effs,
    mixed_r = mixed_r,
    aug_is_residual = aug_is_residual,
    aug_is_numeric = aug_is_numeric,
    aug_error = aug_error,
    df_terms = aug_df_terms,
    den_dfs = den_dfs,
    den_labels = den_labels,
    den_rules = den_rules,
    use_ED = use_ED_active,
    EDall = EDall
  )
  LS_regression_out <- if(isTRUE(add_error)) aug_regression$LS_aug else LS_regression
  LS_SStype_out <- if(isTRUE(add_error)) aug_SStype$LS_aug else LS_SStype
  effect_source <- if(isTRUE(respect_SStype)) "SStype" else "regression"
  LS <- if(isTRUE(respect_SStype)) LS_SStype_out else LS_regression_out
  error <- if(isTRUE(respect_SStype)) aug_SStype$error else aug_regression$error
  LS_aug <- if(isTRUE(respect_SStype)) aug_SStype$LS_aug else aug_regression$LS_aug

  approved <- seq_len(n_terms)
  names(approved) <- effs
  approvedAB <- approved
  approvedComb <- as.list(approved)
  names(approvedComb) <- names(approvedAB)
  eff_combined <- setNames(rep(FALSE, n_terms), effs)

  if(is.logical(REML)){
    if(exists(".ML_variance_partition_single", mode = "function", inherits = TRUE)){
      ssq <- stats::setNames(rep(0, length(effs) + 1), c(effs, "Residuals"))
      for(i in seq_along(models))
        ssq <- ssq + .ML_variance_partition_single(models[[i]], SStype)$Variance * N
    } else if("HDANOVA" %in% loadedNamespaces()){
      vp_fun_post <- get(".ML_variance_partition_single", envir = asNamespace("HDANOVA"))
      ssq <- stats::setNames(rep(0, length(effs) + 1), c(effs, "Residuals"))
      for(i in seq_along(models))
        ssq <- ssq + vp_fun_post(models[[i]], SStype)$Variance * N
    } else {
      warning("Variance-partition helper not found; using QR-based SSQ fallback for REML branch.")
      ssq_terms <- rowSums(SS_matrix)
      ssq <- c(stats::setNames(ssq_terms, effs), Residuals = sum(SSE_full))
    }
  } else {
    ssq_terms <- rowSums(SS_matrix)
    ssq <- c(stats::setNames(ssq_terms, effs), Residuals = sum(SSE_full))
  }
  ssqY <- sum((Y - rep(colMeans(Y), each = N))^2)

  dfNum <- c(df_terms, Residuals = df_res)
  dfDenom <- c(stats::setNames(den_dfs, effs), Residuals = 0)
  denoms <- c(vapply(den_labels, function(lbl){
    if(identical(lbl, "Residuals"))
      return(length(effs) + 1)
    if(lbl %in% effs)
      return(match(lbl, effs))
    NA_real_
  }, numeric(1)), NA_real_)
  names(denoms) <- names(dfDenom)

  model <- model.frame(mod)
  model[[1]] <- Yorig

  obj <- list(LS = LS,
              LS_regression = LS_regression_out,
              LS_SStype = LS_SStype_out,
              effects = effects,
              coefficients = coefs,
              Y = Yorig,
              Y_fit = Y,
              X = M,
              residuals = residuals,
              error = error,
              error_regression = aug_regression$error,
              error_SStype = aug_SStype$error,
              eff_combined = eff_combined,
              SStype = SStype,
              contrasts = contrasts,
              unrestricted = unrestricted,
              ssq = ssq,
              ssqY = ssqY,
              explvar = ssq / ssqY,
              models = if(is.null(models)) list(mod) else models,
              anovas = anovas,
              model.frame = modFra,
              call = match.call(),
              fit.type = fit.type,
              add_error = add_error,
              dfNum = dfNum,
              dfDenom = dfDenom,
              denoms = denoms,
              model = model,
              permute = NULL,
              more = list(approved = approved,
                          approvedAB = approvedAB,
                          effs = effs,
                          assign = assign,
                          approvedComb = approvedComb,
                          N = N,
                          LS_aug = LS_aug,
                          LS_aug_regression = aug_regression$LS_aug,
                          LS_aug_SStype = aug_SStype$LS_aug,
                          REML = REML,
                          lme4 = FALSE,
                          maxDir = maxDir,
                          p = p,
                          pca.in = pca.in,
                          Y_fit = Y,
                          weights = weights_vec,
                          respect_SStype = respect_SStype,
                          effect_source = effect_source,
                          version = 2))
  if(is.logical(REML) && exists(".ML_variance_partition_single", mode = "function", inherits = TRUE)){
    current_options <- options(old_options)
    ssq_final <- stats::setNames(rep(0, length(obj$more$effs) + 1), c(obj$more$effs, "Residuals"))
    for(i in seq_along(obj$models))
      ssq_final <- ssq_final + .ML_variance_partition_single(obj$models[[i]], SStype)$Variance * obj$more$N
    options(current_options)
    obj$ssq <- ssq_final
    obj$explvar <- obj$ssq / obj$ssqY
  }
  if(!is.null(EDall))
    obj$ED <- EDall
  if(pca.in != 0)
    obj$Ypca <- list(pca = pca, ncomp = pca.in)
  class(obj) <- c("hdanova2", "hdanova", "list")
  obj
}


.hda_build_augmented_effect_matrices <- function(LS_base,
                                                 residuals,
                                                 effs,
                                                 mixed_r,
                                                 aug_is_residual,
                                                 aug_is_numeric,
                                                 aug_error,
                                                 df_terms,
                                                 den_dfs,
                                                 den_labels,
                                                 den_rules = NULL,
                                                 use_ED = FALSE,
                                                 EDall = NULL){
  n_terms <- length(effs)
  LS_aug <- vector("list", n_terms)
  error <- vector("list", n_terms)
  names(LS_aug) <- names(error) <- effs
  force_residual <- aug_is_residual || !mixed_r

  for(i in seq_len(n_terms)){
    eff <- effs[i]
    den_idx <- NA_real_
    den_indices <- NA_real_
    den_coefs <- 1
    den_residual_only <- FALSE
    if(force_residual){
      err_component <- residuals
    } else {
      if(!is.null(den_rules) && length(den_rules) >= i){
        rule_i <- den_rules[[i]]
        if(!is.null(rule_i$indices) && length(rule_i$indices) > 0)
          den_indices <- as.numeric(rule_i$indices)
        if(!is.null(rule_i$coefficients) && length(rule_i$coefficients) == length(den_indices))
          den_coefs <- as.numeric(rule_i$coefficients)
      }
      if(all(is.na(den_indices))){
        if(!is.null(den_labels) && length(den_labels) >= i && den_labels[i] %in% effs)
          den_idx <- match(den_labels[i], effs)
        if(is.na(den_idx) && !is.null(den_labels) && length(den_labels) >= i && identical(den_labels[i], "Residuals"))
          den_idx <- length(effs) + 1
        den_indices <- den_idx
        den_coefs <- 1
      }
      den_indices <- den_indices[is.finite(den_indices)]
      if(length(den_indices) == 0){
        den_indices <- length(effs) + 1
        den_coefs <- 1
      }
      den_residual_only <- all(den_indices == (length(effs) + 1))

      if(aug_is_numeric){
        df_num_i <- as.numeric(df_terms[i])
        df_den_i <- as.numeric(den_dfs[i])
        use_effective_dimensions <- isTRUE(use_ED) && !is.null(EDall)
        if(use_effective_dimensions){
          den_names <- ifelse(den_indices == (length(effs) + 1), "Residuals", effs[den_indices])
          if(length(den_names) == 1){
            den_ed <- EDall[den_names, ]
          } else {
            den_ed <- EDall[den_names, , drop = FALSE]
          }
          if(any(is.nan(EDall[eff, ])) || any(is.nan(den_ed))){
            warning("Indeterminate effective dimensions, using degrees of freedom instead")
            use_effective_dimensions <- FALSE
          }
        }
        if(use_effective_dimensions){
          scale_i <- sqrt(EDall[eff, ] / den_ed *
            stats::qf(1 - aug_error, EDall[eff, ], den_ed))
        } else {
          scale_i <- sqrt(df_num_i / df_den_i * stats::qf(1 - aug_error, df_num_i, df_den_i))
        }
      } else {
        scale_i <- 1
      }
      scale_i <- matrix(scale_i, nrow = nrow(residuals), ncol = ncol(residuals), byrow = TRUE)

      err_component <- 0 * residuals
      for(j in seq_along(den_indices)){
        idx_j <- den_indices[j]
        coef_j <- if(length(den_coefs) >= j) den_coefs[j] else 1
        if(idx_j == (length(effs) + 1)){
          err_component <- err_component + coef_j * residuals * scale_i
        } else if(idx_j >= 1 && idx_j <= length(effs)){
          err_component <- err_component + coef_j * LS_base[[effs[idx_j]]] * scale_i
        }
      }
      den_idx <- if(length(den_indices) == 1) den_indices else NA_real_
    }

    if(!force_residual && !is.na(den_idx) && den_idx <= length(effs)){
      error[[eff]] <- LS_base[[eff]] + err_component
      LS_aug[[eff]] <- LS_base[[eff]] + error[[eff]]
    } else {
      LS_aug[[eff]] <- LS_base[[eff]] + err_component
      error[[eff]] <- LS_aug[[eff]]
      if(!force_residual && isTRUE(den_residual_only))
        attr(error[[eff]], "term") <- "residuals"
    }
  }

  list(LS_aug = LS_aug, error = error)
}


.hda_compute_effective_dimensions <- function(models, Y, X, Z, df_template){
  varcor_random <- lapply(models, lme4::VarCorr)
  varcor_resid <- vapply(models, function(mod) stats::sigma(mod)^2, numeric(1))
  randEffList <- lme4::getME(models[[1]], "flist")
  if(length(randEffList) == 0)
    return(NULL)

  # TODO: This ED calculation was copied from another package. Keep parity for
  # now, but revisit the derivation in a dedicated ED review.
  ED <- vector("list", ncol(Y))
  N <- nrow(Y)
  q <- ncol(Z)

  for(i in seq_len(ncol(Y))){
    wrong.ED <- FALSE
    G <- diag(q) * 0
    r0 <- 0
    for(j in seq_along(randEffList)){
      # TODO: Revisit broader random-effect structures later. This currently
      # assumes full coding of random effects.
      nlevR <- nlevels(randEffList[[j]])
      G[r0 + (1:nlevR), r0 + (1:nlevR)] <- diag(nlevR) *
        max(varcor_random[[i]][[names(randEffList)[j]]][1], 10^-12)
      r0 <- r0 + nlevR
    }
    Rinv <- diag(N) / varcor_resid[i]

    XRX <- t(X) %*% Rinv %*% X
    ZRX <- t(Z) %*% Rinv %*% X
    ZRZ <- t(Z) %*% Rinv %*% Z
    Ginv <- NULL
    Ginv_try <- try(solve(G), silent = TRUE)
    if(inherits(Ginv_try, "try-error")){
      Ginv <- G * Inf
      wrong.ED <- TRUE
      warning("Random effects covariance matrix is singular, using Inf instead")
    } else {
      Ginv <- Ginv_try
    }
    ZRZG <- ZRZ + Ginv
    Q <- solve(rbind(cbind(XRX, t(ZRX)), cbind(ZRX, ZRZG)))
    K <- Q %*% cbind(rbind(XRX, ZRX), rbind(t(ZRX), ZRZ))
    K_diag <- diag(K)[-seq_len(ncol(X))]
    ED_rand <- numeric(length(randEffList))
    r0 <- 0
    for(j in seq_along(randEffList)){
      nlevR <- nlevels(randEffList[[j]])
      ED_rand[j] <- sum(K_diag[r0 + (1:nlevR)])
      r0 <- r0 + nlevR
    }
    if(wrong.ED)
      ED_rand <- rep(NaN, length(ED_rand))
    ED[[i]] <- ED_rand
  }

  EDm <- matrix(unlist(ED), nrow = length(randEffList))
  EDm <- pmax(EDm, 1)
  rownames(EDm) <- names(randEffList)

  EDall <- matrix(df_template, nrow = length(df_template), ncol = ncol(EDm))
  rownames(EDall) <- names(df_template)
  colnames(EDall) <- colnames(Y)
  EDall["Residuals", ] <- colSums(EDall) + 1 - colSums(EDm)
  EDall[rownames(EDm), ] <- EDm
  EDall
}


.calc_SS_qr_loop <- function(X, Y, assign, term_labels, SStype, weights = NULL){
  n_terms <- length(term_labels)
  p <- ncol(Y)
  SS_matrix <- matrix(0, nrow = n_terms, ncol = p,
                      dimnames = list(term_labels, colnames(Y)))
  effects <- setNames(vector("list", n_terms), term_labels)
  df_terms <- vapply(seq_len(n_terms), function(i) sum(assign == i), integer(1))
  names(df_terms) <- term_labels
  sqrt_weights <- if(is.null(weights)) rep(1, nrow(X)) else sqrt(weights)

  full_fit <- .weighted_qr_fit(X, Y, sqrt_weights)
  QR_full <- full_fit$QR
  Y_hat_full <- full_fit$fitted

  for(i in seq_len(n_terms)){
    if(SStype == 1){
      Y_hat_curr <- .weighted_qr_fit(X[, which(assign <= i), drop = FALSE], Y, sqrt_weights)$fitted
      Y_hat_prev <- .weighted_qr_fit(X[, which(assign < i), drop = FALSE], Y, sqrt_weights)$fitted
    } else if(SStype == 2){
      is_higher <- grepl(paste0("^", term_labels[i], ":"), term_labels) |
        grepl(paste0(":", term_labels[i], "$"), term_labels)
      Y_hat_curr <- .weighted_qr_fit(X[, which(!(assign %in% which(is_higher))), drop = FALSE], Y, sqrt_weights)$fitted
      Y_hat_prev <- .weighted_qr_fit(X[, which(!(assign %in% c(i, which(is_higher)))), drop = FALSE], Y, sqrt_weights)$fitted
    } else {
      Y_hat_curr <- Y_hat_full
      Y_hat_prev <- .weighted_qr_fit(X[, which(assign != i), drop = FALSE], Y, sqrt_weights)$fitted
    }
    effects[[term_labels[i]]] <- Y_hat_curr - Y_hat_prev
    SS_matrix[i, ] <- colSums((effects[[term_labels[i]]] * sqrt_weights)^2)
  }

  list(SS_matrix = SS_matrix,
       effects = effects,
       df_terms = df_terms,
       QR_full = QR_full,
       Y_hat_full = Y_hat_full)
}


.weighted_qr_fit <- function(X, Y, sqrt_weights){
  X_weighted <- X * sqrt_weights
  Y_weighted <- Y * sqrt_weights
  QR <- qr(X_weighted)
  coefs <- qr.coef(QR, Y_weighted)
  estimable <- !is.na(rowSums(coefs))
  fitted <- matrix(0, nrow = nrow(X), ncol = ncol(Y), dimnames = dimnames(Y))
  if(any(estimable))
    fitted <- X[, estimable, drop = FALSE] %*% coefs[estimable, , drop = FALSE]
  list(QR = QR, coefficients = coefs, fitted = fitted)
}
