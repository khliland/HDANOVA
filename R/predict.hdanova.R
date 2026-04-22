#' Predict for HDANOVA Objects
#'
#' Reconstructs an HDANOVA-style object on \code{newdata} without refitting by
#' reusing stored coefficients and projection objects from the fitted model.
#' This implementation supports fixed-effects, mixed MoM (\code{r()} with
#' \code{REML = NULL}), and REML/ML mixed workflows (\code{r()} with
#' \code{REML = TRUE} or \code{FALSE}).
#'
#' @param object A fitted \code{hdanova} object.
#' @param newdata A data frame containing the variables used in the original model formula.
#' @param ... Reserved for generic compatibility; runtime overrides are not supported.
#'
#' @return An \code{hdanova}-family object computed on \code{newdata}.
#'
#' @seealso ASCA-family prediction wrappers: \code{\link{predict_asca_family}}.
#' Model constructors: \code{\link{hdanova}}, \code{\link{asca}},
#' \code{\link{apca}}, \code{\link{apls}}, \code{\link{limmpca}},
#' \code{\link{msca}} and \code{\link{pcanova}}.
#'
#' @examples
#' data(candies)
#' # Train/test split (every third sample to test)
#' test_idx  <- seq(3, nrow(candies), by = 3)
#' train_idx <- setdiff(1:nrow(candies),test_idx)
#' candies_train <- candies[train_idx, ]
#' candies_test  <- candies[test_idx, ]
#'
#' # Fixed-effects model prediction
#' mod <- hdanova(assessment ~ candy + assessor, data = candies_train)
#' pred <- predict(mod, newdata = candies_test)
#'
#' var_idx <- seq_len(ncol(mod$LS$candy))
#' old.par <- par(mfrow = c(1,2), mar = c(4,4,2,1), mgp = c(2,0.7,0))
#' image(x = var_idx, y = seq_along(train_idx), z = t(mod$LS$candy),
#'       xaxt = "n", yaxt = "n", main = "Original candy LS",
#'       xlab = "Variable index", ylab = "Train sample index")
#' axis(1, at = var_idx, labels = var_idx)
#' axis(2, at = seq_along(train_idx), labels = train_idx)
#' image(x = var_idx, y = seq_along(test_idx), z = t(pred$LS$candy),
#'       xaxt = "n", yaxt = "n", main = "Predicted candy LS",
#'       xlab = "Variable index", ylab = "Test sample index")
#' axis(1, at = var_idx, labels = var_idx)
#' axis(2, at = seq_along(test_idx), labels = test_idx)
#' par(old.par)
#'
#' # Mixed MoM model prediction (r() with REML = NULL)
#' mod_mom <- hdanova(assessment ~ candy + r(assessor), data = candies_train)
#' pred_mom <- predict(mod_mom, newdata = candies_test)
#' cat("Mixed MoM model prediction successful.\n")
#' cat("SSQ names:", paste(names(pred_mom$ssq), collapse = "|"), "\n")
#' cat("dfDenom:", paste(pred_mom$dfDenom, collapse = "|"), "\n")
#'
#' # REML mixed model prediction (r() with REML = TRUE)
#' mod_reml <- hdanova(assessment ~ candy + r(assessor), data = candies_train, REML = TRUE)
#' pred_reml <- predict(mod_reml, newdata = candies_test)
#' cat("REML mixed model prediction successful.\n")
#' cat("SSQ names:", paste(names(pred_reml$ssq), collapse = "|"), "\n")
#' cat("dfDenom:", paste(pred_reml$dfDenom, collapse = "|"), "\n")
#'
#' @method predict hdanova
#' @export
predict.hdanova <- function(object, newdata, ...){
  if(missing(newdata))
    stop("'newdata' must be provided.")
  if(!is.data.frame(newdata))
    stop("'newdata' must be a data.frame.")

  dots <- list(...)
  if(length(dots) > 0)
    stop("predict.hdanova() does not support runtime overrides; use training-time settings.")

  call_obj <- object[["call"]]
  if(is.null(call_obj) || !inherits(call_obj, "call"))
    stop("The supplied object does not contain a reproducible model call.")
  if(is.null(call_obj[["formula"]]))
    stop("Could not extract the original formula from 'object$call'.")

  formula_txt <- paste(deparse(call_obj[["formula"]]), collapse = " ")
  has_random_terms <- grepl("r(", formula_txt, fixed = TRUE)
  mixed_mom <- isTRUE(has_random_terms) && is.null(object[["more"]][["REML"]])
  mixed_reml <- isTRUE(has_random_terms) && is.logical(object[["more"]][["REML"]])
  mixed_r <- isTRUE(mixed_mom) || isTRUE(mixed_reml)

  clean_formula_from_r_terms <- function(formula_obj, effs_fixed){
    resp_node <- formula_obj[[2]]
    rhs_terms <- paste(effs_fixed, collapse = " + ")
    stats::as.formula(paste(paste(resp_node, "~", sep=" "), rhs_terms))
  }

  form_fixed <- if(isTRUE(mixed_r)){
    effs_fixed <- object[["more"]][["effs_fixed"]]
    if(is.null(effs_fixed) || length(effs_fixed) == 0){
      clean_formula_from_r_terms(call_obj[["formula"]], object$more$effs)
    } else {
      clean_formula_from_r_terms(call_obj[["formula"]], effs_fixed)
    }
  } else {
    stats::formula(object$models[[1]])
  }

  train_models <- object[["models"]]
  if(!is.list(train_models) || length(train_models) == 0)
    stop("No fitted model information found in 'object$models'.")
  train_mod <- train_models[[1]]

  xlevels <- NULL
  if(inherits(train_mod, "lmerMod")){
    if(!is.null(try(train_mod@frame, silent = TRUE))){
      frame <- train_mod@frame
      xlevels <- list()
      for(vn in names(frame)){
        if(is.factor(frame[[vn]])){
          xlevels[[vn]] <- levels(frame[[vn]])
        }
      }
    }
  } else {
    xlevels <- train_mod[["xlevels"]]
  }
  
  if(is.list(xlevels) && length(xlevels) > 0){
    for(vn in intersect(names(xlevels), names(newdata))){
      vals <- as.character(newdata[[vn]])
      unseen <- setdiff(unique(vals[!is.na(vals)]), xlevels[[vn]])
      if(length(unseen) > 0){
        stop(paste0("newdata contains unseen level(s) for '", vn, "': ",
                    paste(unseen, collapse = ", "), "."))
      }
      newdata[[vn]] <- factor(vals, levels = xlevels[[vn]])
    }
  }

  mf_data <- stats::model.frame(form_fixed, data = newdata, na.action = stats::na.pass)
  Y <- stats::model.response(mf_data)
  if(inherits(Y, "AsIs"))
    Y <- unclass(Y)
  Yorig <- Y <- as.matrix(Y)
  if(is.null(ncol(Y)))
    stop("Response in 'newdata' must be a matrix-like block.")

  pca_in <- object[["more"]][["pca.in"]]
  pls_in <- object[["more"]][["pls.in"]]
  pca_active <- !(is.logical(pca_in) && !pca_in) && !identical(pca_in, 0)
  pls_active <- !(is.logical(pls_in) && !pls_in) && !identical(pls_in, 0)
  if(pca_active && pls_active)
    stop("Object has both active PCA and PLS compression; cannot predict reliably.")

  if(pca_active){
    if(is.null(object[["Ypca"]]) || is.null(object[["Ypca"]][["pca"]]))
      stop("Object indicates pca.in but does not contain stored PCA model in 'object$Ypca'.")
    pca_obj <- object[["Ypca"]][["pca"]]
    ncomp <- object[["Ypca"]][["ncomp"]]
    loads <- as.matrix(pca_obj[["loadings"]])[, seq_len(ncomp), drop = FALSE]
    means <- as.numeric(pca_obj[["Xmeans"]])
    if(length(means) != ncol(Y))
      stop("Stored PCA mean vector length does not match response width in 'newdata'.")
    Y <- (Y - rep(means, each = nrow(Y))) %*% loads
  } else if(pls_active){
    if(is.null(object[["Ypls"]]) || is.null(object[["Ypls"]][["pls"]]))
      stop("Object indicates pls.in but does not contain stored PLS model in 'object$Ypls'.")
    pls_obj <- object[["Ypls"]][["pls"]]
    ncomp <- object[["Ypls"]][["ncomp"]]
    proj <- pls_obj[["projection"]]
    means <- as.numeric(pls_obj[["Xmeans"]])
    if(is.null(proj) || length(means) != ncol(Y))
      stop("Stored PLS projection is incompatible with response width in 'newdata'.")
    Y <- (Y - rep(means, each = nrow(Y))) %*% proj[, seq_len(ncomp), drop = FALSE]
  }

  tt <- stats::delete.response(stats::terms(train_mod))
  
  contrasts_arg <- NULL
  if(inherits(train_mod, "lmerMod")){
    mm_train_full <- stats::model.matrix(train_mod)
    contrasts_arg <- attr(mm_train_full, "contrasts")
  } else {
    contrasts_arg <- train_mod[["contrasts"]]
  }
  
  M_raw <- stats::model.matrix(tt, data = mf_data, contrasts.arg = contrasts_arg)
  assign_raw <- attr(M_raw, "assign")
  effs_base <- attr(tt, "term.labels")
  n_terms <- length(effs_base)
  if(n_terms == 0)
    stop("No fixed effects found in model terms.")

  coefs <- as.matrix(object[["coefficients"]])
  if(ncol(coefs) != ncol(Y))
    stop("Stored coefficient dimensions are incompatible with projected response dimensions.")

  if(isTRUE(mixed_reml)){
    X_fixed_stored <- object[["more"]][["X_fixed"]]
    if(!is.null(X_fixed_stored)){
      n_fixed_cols <- ncol(X_fixed_stored)
      coefs <- coefs[seq_len(n_fixed_cols), , drop = FALSE]
    }
  }

  coef_names <- rownames(coefs)
  mm_names <- colnames(M_raw)

  coef_to_mm_name <- function(coef_name, xlevels){
    if(identical(coef_name, "(Intercept)"))
      return(coef_name)
    pieces <- strsplit(coef_name, ":", fixed = TRUE)[[1]]
    out <- vapply(pieces, function(pc){
      mt <- regexec("^(.+)\\((.*)\\)$", pc)
      rg <- regmatches(pc, mt)[[1]]
      if(length(rg) == 3){
        var <- rg[2]
        lev <- rg[3]
        if(!is.null(xlevels[[var]])){
          idx <- match(lev, xlevels[[var]])
          if(!is.na(idx))
            return(paste0(var, idx))
        }
      }
      pc
    }, character(1))
    paste(out, collapse = ":")
  }

  map <- match(coef_names, mm_names)
  if(anyNA(map)){
    mapped_names <- vapply(coef_names, coef_to_mm_name, character(1), xlevels = xlevels)
    map <- match(mapped_names, mm_names)
  }
  if(anyNA(map) || any(duplicated(map)))
    stop("Coefficient/design column mismatch between object and newdata model matrix.")

  M <- M_raw[, map, drop = FALSE]
  colnames(M) <- coef_names
  assign <- assign_raw[map]
  coefs_aligned <- coefs

  LS_regression <- vector("list", n_terms)
  names(LS_regression) <- effs_base
  for(i in seq_len(n_terms)){
    idx <- which(assign == i)
    Li <- matrix(0, nrow = nrow(Y), ncol = ncol(Y), dimnames = list(rownames(Y), colnames(Y)))
    if(length(idx) > 0){
      ci <- coefs_aligned[idx, , drop = FALSE]
      estimable <- !is.na(rowSums(ci))
      if(any(estimable)){
        Li <- M[, idx, drop = FALSE][, estimable, drop = FALSE] %*%
          ci[estimable, , drop = FALSE]
      }
    }
    LS_regression[[effs_base[i]]] <- Li
  }

  SStype <- object[["SStype"]]
  qr_res <- .calc_SS_qr_loop(M, Y, assign, effs_base, SStype, weights = NULL)
  LS_SStype <- qr_res$effects

  estimable_all <- !is.na(rowSums(coefs_aligned))
  Y_hat <- matrix(0, nrow = nrow(Y), ncol = ncol(Y), dimnames = dimnames(Y))
  if(any(estimable_all))
    Y_hat <- M[, estimable_all, drop = FALSE] %*% coefs_aligned[estimable_all, , drop = FALSE]
  residuals <- Y - Y_hat

  SS_matrix <- qr_res$SS_matrix
  df_terms <- qr_res$df_terms
  df_res <- nrow(Y) - qr_res$QR_full$rank
  SSE_full <- colSums(residuals^2)
  MS_matrix <- sweep(SS_matrix, 1, df_terms, "/")
  residual_df_valid <- df_res > 0
  MS_res <- if(residual_df_valid) SSE_full / df_res else rep(NA_real_, ncol(Y))

  den_dfs <- stats::setNames(rep(if(residual_df_valid) df_res else NA_real_, n_terms), effs_base)
  den_labels <- stats::setNames(rep("Residuals", n_terms), effs_base)
  MS_denom <- matrix(NA_real_, nrow = n_terms, ncol = ncol(MS_matrix),
                     dimnames = dimnames(MS_matrix))
  if(residual_df_valid)
    MS_denom[] <- rep(MS_res, each = n_terms)

  train_df_denom <- object[["dfDenom"]]
  if(!is.null(train_df_denom)){
    shared_den_df <- intersect(effs_base, names(train_df_denom))
    den_dfs[shared_den_df] <- as.numeric(train_df_denom[shared_den_df])
  }

  train_denoms <- object[["denoms"]]
  if(!is.null(train_denoms)){
    shared_den <- intersect(effs_base, names(train_denoms))
    for(eff in shared_den){
      den_idx <- suppressWarnings(as.integer(train_denoms[[eff]]))
      if(is.na(den_idx) || den_idx == (n_terms + 1)){
        den_labels[[eff]] <- "Residuals"
      } else if(den_idx >= 1 && den_idx <= n_terms){
        den_labels[[eff]] <- effs_base[den_idx]
        MS_denom[eff, ] <- MS_matrix[den_idx, ]
      }
    }
  }

  F_matrix <- P_matrix <- matrix(NA_real_, nrow = n_terms, ncol = ncol(Y),
                                 dimnames = list(effs_base, colnames(Y)))
  for(i in seq_len(n_terms)){
    valid_cols <- is.finite(MS_matrix[i, ]) & is.finite(MS_denom[i, ]) &
      !isTRUE(all.equal(den_dfs[i], 0)) & is.finite(den_dfs[i])
    if(any(valid_cols)){
      F_matrix[i, valid_cols] <- MS_matrix[i, valid_cols] / MS_denom[i, valid_cols]
      P_matrix[i, valid_cols] <- stats::pf(F_matrix[i, valid_cols], df_terms[i], den_dfs[i], lower.tail = FALSE)
    }
  }

  anovas <- vector("list", ncol(Y))
  names(anovas) <- colnames(Y)
  for(i in seq_len(ncol(Y))){
    ano <- data.frame("Df" = c(df_terms, df_res),
                      "Sum Sq" = c(SS_matrix[, i], SSE_full[i]),
                      "Mean Sq" = c(MS_matrix[, i], MS_res[i]),
                      "F value" = c(F_matrix[, i], NA_real_),
                      "Pr(>F)" = c(P_matrix[, i], NA_real_),
                      "Error Term" = c(den_labels, NA_character_),
                      check.names = FALSE)
    rownames(ano) <- c(effs_base, "Residuals")
    anovas[[i]] <- ano
  }

  modFra <- HDANOVA::extended.model.frame(form_fixed, mf_data)
  effects <- setNames(lapply(effs_base, function(eff) modFra[[eff]]), effs_base)

  respect_SStype <- object[["more"]][["respect_SStype"]]
  add_error <- object[["add_error"]]

  aug_error <- "denominator"
  if(!is.null(call_obj[["aug_error"]]))
    aug_error <- eval(call_obj[["aug_error"]], envir = parent.frame())
  aug_is_residual <- is.character(aug_error) && length(aug_error) == 1 &&
    aug_error %in% c("residual", "residuals")
  aug_is_numeric <- is.numeric(aug_error) && length(aug_error) == 1 && is.finite(aug_error)

  aug_regression <- .hda_build_augmented_effect_matrices(
    LS_base = LS_regression,
    residuals = residuals,
    effs = effs_base,
    mixed_r = mixed_r,
    aug_is_residual = aug_is_residual,
    aug_is_numeric = aug_is_numeric,
    aug_error = aug_error,
    df_terms = df_terms,
    den_dfs = den_dfs,
    den_labels = den_labels,
    den_rules = NULL,
    use_ED = FALSE,
    EDall = NULL
  )
  aug_SStype <- .hda_build_augmented_effect_matrices(
    LS_base = LS_SStype,
    residuals = residuals,
    effs = effs_base,
    mixed_r = mixed_r,
    aug_is_residual = aug_is_residual,
    aug_is_numeric = aug_is_numeric,
    aug_error = aug_error,
    df_terms = df_terms,
    den_dfs = den_dfs,
    den_labels = den_labels,
    den_rules = NULL,
    use_ED = FALSE,
    EDall = NULL
  )

  combine_effect_list <- function(base_list){
    out <- base_list
    comb_flags <- object[["eff_combined"]]
    comb_map <- object[["more"]][["approvedComb"]]
    if(is.null(comb_flags) || is.null(comb_map))
      return(out)
    comb_names <- names(comb_flags)[isTRUE(comb_flags)]
    for(cn in comb_names){
      members_idx <- comb_map[[cn]]
      if(is.null(members_idx))
        next
      members <- names(members_idx)
      if(is.null(members) || any(members == ""))
        members <- effs_base[as.integer(members_idx)]
      members <- members[members %in% names(base_list)]
      if(length(members) == 0)
        next
      out[[cn]] <- Reduce("+", lapply(members, function(mn) base_list[[mn]]))
    }
    out
  }

  LS_regression_all <- combine_effect_list(LS_regression)
  LS_SStype_all <- combine_effect_list(LS_SStype)
  LS_aug_reg_all <- combine_effect_list(aug_regression$LS_aug)
  LS_aug_ss_all <- combine_effect_list(aug_SStype$LS_aug)
  err_reg_all <- combine_effect_list(aug_regression$error)
  err_ss_all <- combine_effect_list(aug_SStype$error)

  LS_regression_out <- if(isTRUE(add_error)) LS_aug_reg_all else LS_regression_all
  LS_SStype_out <- if(isTRUE(add_error)) LS_aug_ss_all else LS_SStype_all
  LS <- if(isTRUE(respect_SStype)) LS_SStype_out else LS_regression_out
  error <- if(isTRUE(respect_SStype)) err_ss_all else err_reg_all
  LS_aug <- if(isTRUE(respect_SStype)) LS_aug_ss_all else LS_aug_reg_all

  ssq_base <- c(stats::setNames(rowSums(SS_matrix), effs_base), Residuals = sum(SSE_full))
  approved <- object[["more"]][["approved"]]
  if(is.null(approved)){
    approved <- seq_len(length(effs_base))
    names(approved) <- effs_base
  }
  approved_names <- names(approved)
  if(is.null(approved_names) || length(approved_names) == 0)
    approved_names <- effs_base
  comb_map <- object[["more"]][["approvedComb"]]
  eff_combined <- object[["eff_combined"]]
  if(is.null(eff_combined))
    eff_combined <- stats::setNames(rep(FALSE, length(approved_names)), approved_names)

  ssq <- stats::setNames(rep(0, length(approved_names) + 1), c(approved_names, "Residuals"))
  for(nm in approved_names){
    if(isTRUE(eff_combined[[nm]]) && !is.null(comb_map[[nm]])){
      members <- names(comb_map[[nm]])
      if(is.null(members) || any(members == ""))
        members <- effs_base[as.integer(comb_map[[nm]])]
      members <- intersect(members, names(ssq_base))
      ssq[[nm]] <- sum(ssq_base[members])
    } else if(nm %in% names(ssq_base)){
      ssq[[nm]] <- ssq_base[[nm]]
    }
  }
  ssq[["Residuals"]] <- ssq_base[["Residuals"]]

  dfNum_eff <- stats::setNames(rep(NA_real_, length(approved_names)), approved_names)
  for(nm in approved_names){
    if(isTRUE(eff_combined[[nm]]) && !is.null(comb_map[[nm]])){
      members <- names(comb_map[[nm]])
      if(is.null(members) || any(members == ""))
        members <- effs_base[as.integer(comb_map[[nm]])]
      members <- intersect(members, names(df_terms))
      dfNum_eff[[nm]] <- sum(df_terms[members])
    } else if(nm %in% names(df_terms)){
      dfNum_eff[[nm]] <- df_terms[[nm]]
    }
  }
  dfNum <- c(dfNum_eff, Residuals = df_res)
  dfDenom <- c(stats::setNames(rep(df_res, length(approved_names)), approved_names), Residuals = 0)
  denoms <- c(stats::setNames(rep(length(approved_names) + 1, length(approved_names)), approved_names), NA_real_)
  names(denoms) <- names(dfDenom)

  train_df_denom <- object[["dfDenom"]]
  if(!is.null(train_df_denom)){
    shared_out_df <- intersect(approved_names, names(train_df_denom))
    dfDenom[shared_out_df] <- as.numeric(train_df_denom[shared_out_df])
  }

  train_denoms <- object[["denoms"]]
  if(!is.null(train_denoms)){
    shared_out_den <- intersect(approved_names, names(train_denoms))
    denoms[shared_out_den] <- as.numeric(train_denoms[shared_out_den])
  }

  ssqY <- sum((Y - rep(colMeans(Y), each = nrow(Y)))^2)
  explvar <- ssq / ssqY

  model <- mf_data
  model[[1]] <- Yorig

  out <- object
  out$LS <- LS
  out$LS_regression <- LS_regression_out
  out$LS_SStype <- LS_SStype_out
  out$effects <- effects
  out$Y <- Yorig
  out$Y_fit <- Y
  out$X <- M
  out$residuals <- residuals
  out$error <- error
  out$error_regression <- err_reg_all
  out$error_SStype <- err_ss_all
  out$eff_combined <- eff_combined
  out$SStype <- SStype
  out$ssq <- ssq
  out$ssqY <- ssqY
  out$explvar <- explvar
  out$anovas <- anovas
  out$model.frame <- modFra
  out$call <- match.call()
  out$add_error <- add_error
  out$dfNum <- dfNum
  out$dfDenom <- dfDenom
  out$denoms <- denoms
  out$model <- model
  out$permute <- NULL

  out$more$assign <- assign
  out$more$assign_fixed <- assign
  out$more$X_fixed <- M
  out$more$effs_fixed <- effs_base
  out$more$N <- nrow(Y)
  out$more$LS_aug <- LS_aug
  out$more$LS_aug_regression <- LS_aug_reg_all
  out$more$LS_aug_SStype <- LS_aug_ss_all
  out$more$maxDir <- vapply(out$more$effs, function(eff){
    if(!eff %in% names(out$LS))
      return(0)
    min(pracma::Rank(out$LS[[eff]]), ncol(Y))
  }, numeric(1))
  out$more$p <- ncol(Y)
  out$more$Y_fit <- Y
  out$more$weights <- NULL
  out$more$respect_SStype <- respect_SStype
  out$more$effect_source <- if(isTRUE(respect_SStype)) "SStype" else "regression"
  out$more$ssq_method <- if(isTRUE(respect_SStype)) "qr_sstype" else "qr_regression"
  out$more$ssq_source <- "qr"

  attr(out, "predicted_from") <- call_obj
  out
}
