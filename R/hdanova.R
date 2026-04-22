#' @name hdanova
#' @aliases hdanova
#' @title High-Dimensional Analysis of Variance
#'
#' @description This function provides a high-dimensional analysis of variance (HDANOVA) method
#' which can be used alone or as part of a larger analysis, e.g., ASCA, APCA, LiMM-PCA, MSCA or PC-ANOVA. It
#' can be called directly or through the convenience functions \code{\link{asca}}, \code{\link{apca}},
#' \code{\link{limmpca}}, \code{\link{msca}} and \code{\link{pcanova}}.
#'
#' @param formula Model formula accepting a single response (block) and predictors. See Details for more information.
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param weights Optional object weights.
#' @param na.action How to handle NAs (no action implemented).
#' @param family Error distributions and link function for Generalized Linear Models.
#' @param scale Scaling of the response matrix. Defaults to \code{FALSE} (no scaling). For alternatives, see Details.
#' @param add_error Add error to LS means, e.g., for APCA.
#' @param aug_error Augment score matrices in backprojection. Default = "denominator"
#' (of F test), "residual" (force error term), numeric value (alpha-value in LiMM-PCA).
#' @param pca.in Compress response before ASCA (number of components).
#' @param pls.in Compress response before ASCA using PLS (number of components).
#' @param contrasts Effect coding: "contr.sum" (default = sum-coding), "contr.weighted" (not for lme4 models), "contr.reference", "contr.treatment".
#' @param unrestricted Use unrestricted ANOVA decomposition (default = FALSE).
#' @param SStype Type of sum-of-squares: "I" = sequential, "II" (default) = last term, obeying marginality,
#' "III" = last term, not obeying marginality.
#' @param respect_SStype Logical; if \code{FALSE} (default), keep regression-based
#' effect matrices in \code{LS}. If \code{TRUE}, expose SS-type-aligned
#' effect matrices in \code{LS} while retaining regression matrices in
#' \code{LS_regression}. This setting also propagates to
#' \code{permutation()} when \code{respect_SStype = NULL} is used there.
#' @param REML Parameter to mixlm: NULL (default) = sum-of-squares, TRUE = REML, FALSE = ML.
#' If supplied without any \code{r()} terms in the formula, it is ignored with a warning and
#' the model is fitted as fixed-effects.
#' @param REML_ssq_method Method for REML mixed-model SSQ decomposition:
#' \code{"exact_refit"} (default), \code{"wald"}, or \code{"ls"}. This is
#' only used when \code{REML} is \code{TRUE} or \code{FALSE} for mixed models
#' with \code{r()} terms.
#' @param equal_baseline Experimental: Set to \code{TRUE} to let interactions, where a main effect is missing,
#' e.g., a nested model, be handled with the same baseline as a cross effect model. If \code{TRUE} the corresponding
#' interactions will be put in quotation marks and included in the \code{model.frame}.
#' @param use_ED Use "effective dimensions" for score rescaling in LiMM-PCA.
#'
#' @details Scaling of the response matrix can be done by setting the \code{scale} parameter. If \code{scale=TRUE},
#' each column is scaled by its standard deviation (autoscaling). A numeric value can be provided to scale
#' the columns by specific quantities. If \code{scale} is a character string, the first element
#' is interpreted as a factor name and the second element is interpreted as a factor level, whose samples
#' the standard deviations are calculated from (reference group scaling).
#'
#' \strong{SStype and respect_SStype (fixed and MoM mixed models):}
#' \code{SStype} controls how each effect's contribution is isolated from the others.
#' Type I (sequential) assigns variance in the order terms appear in the formula;
#' Type II (default) tests each term against all others that do not contain it, respecting
#' marginality; Type III tests each term against the full model, ignoring marginality.
#' For balanced designs the three types give identical results; differences arise in
#' unbalanced or non-orthogonal designs.
#'
#' By default (\code{respect_SStype = FALSE}) the effect matrices in \code{LS} are regression
#' projections — each effect's columns are projected independently from the full model
#' coefficient matrix. This is consistent with the legacy ASCA workflow. When
#' \code{respect_SStype = TRUE} the \code{LS} matrices are re-derived from QR contrasts
#' that match the chosen \code{SStype}, so the Frobenius norm of \code{LS[[eff]]} equals
#' the corresponding ANOVA sum of squares. The regression matrices are still available in
#' \code{LS_regression}. For balanced, fully-crossed designs the two sets are numerically
#' identical; for unbalanced designs they will differ, particularly for Type II and III.
#'
#' \strong{REML SSQ strategies (mixed models with \code{r()} terms):}
#' When \code{REML = TRUE} or \code{REML = FALSE}, the model is fitted with lme4 and the
#' classical ANOVA table is replaced by one of three decomposition strategies controlled by
#' \code{REML_ssq_method}:
#' \itemize{
#'   \item \code{"exact_refit"} (default): For each effect a reduced model is re-fitted
#'     (REML or ML as specified) and the SSQ is the difference in log-likelihoods scaled
#'     to the sum-of-squares metric. Concretely, if \eqn{\ell_f} and \eqn{\ell_r} are the
#'     log-likelihoods of the full and reduced models, the statistic is
#'     \deqn{\mathrm{SSQ} = 2(\ell_f - \ell_r).}
#'     This is the most principled approach but requires
#'     one additional model fit per effect, making it slower for large models.
#'   \item \code{"wald"}: Uses the Wald chi-square statistic from the full model divided
#'     by the residual variance to obtain an approximate SSQ. For a fixed effect with
#'     coefficient vector \eqn{\hat{\boldsymbol{\beta}}_j} and covariance
#'     \eqn{\mathbf{V}_j = \mathrm{Var}(\hat{\boldsymbol{\beta}}_j)},
#'     \deqn{\mathrm{SSQ} \approx \hat{\boldsymbol{\beta}}_j^\top \mathbf{V}_j^{-1} \hat{\boldsymbol{\beta}}_j \cdot \hat{\sigma}^2.}
#'     No additional model fits are
#'     required; fast but less accurate for small samples or near-singular random effects.
#'   \item \code{"ls"}: Projects the REML/ML coefficient estimates through the fixed-effect
#'     design matrix to recover a least-squares-style SSQ,
#'     \deqn{\mathrm{SSQ} = \| \mathbf{X}_j \hat{\boldsymbol{\beta}}_j \|_F^2,}
#'     where \eqn{\mathbf{X}_j} contains the columns of the design matrix for effect \eqn{j}.
#'     Computationally cheap and
#'     numerically stable but ignores the mixed-model covariance structure; best treated
#'     as an approximation.
#' }
#'
#' \strong{Permutation statistics vs. fitted-model SSQ:}
#' Permutation testing (\code{\link{permutation}}) always uses QR-based projection
#' statistics computed on the fixed-effect design matrix, regardless of \code{REML_ssq_method}.
#' When \code{respect_SStype = FALSE} (default) the permutation statistic is a
#' regression-projection norm; when \code{respect_SStype = TRUE} it is an SS-type-aligned
#' QR contrast norm. Neither is identical to the fitted-model REML SSQ values reported in
#' \code{object$ssq}, because REML/ML decompositions account for the random-effect
#' covariance structure whereas permutation statistics do not. The p-values are still
#' valid under their respective null hypotheses; the SSQ values should not be compared
#' across the two sources.
#'
#' @return An \code{hdanova} object containing loadings, scores, explained variances, etc. The object has
#' associated plotting (\code{\link{asca_plots}}) and result (\code{\link{asca_results}}) functions.
#'
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic HDANOVA model with two factors
#' mod <- hdanova(assessment ~ candy + assessor, data=candies)
#' summary(mod)
#'
#' @importFrom lme4 lmer glmer getME ranef VarCorr fixef
#' @importFrom progress progress_bar
#' @importFrom grDevices adjustcolor palette
#' @importFrom graphics abline axis box hist legend lines points
#' @importFrom stats anova coefficients contr.sum contr.treatment contrasts<- fitted formula getCall logLik model.frame model.matrix model.response qf reformulate rnorm setNames sigma terms update var
#' @importFrom pracma Rank
#' @importFrom MASS ginv
#' @export
hdanova <- function(formula, data, subset, weights, na.action, family,
                    scale = FALSE,
                    add_error = FALSE, # TRUE => APCA
                    aug_error = "denominator", # "residual" => Mixed, alpha-value => LiMM-PCA
                    pca.in = FALSE,
                    pls.in = FALSE,
                    contrasts = "contr.sum",
                    unrestricted = FALSE,
                    SStype = "II",
                    respect_SStype = FALSE,
                    REML = NULL,
                    REML_ssq_method = c("exact_refit", "wald", "ls"),
                    equal_baseline = FALSE,
                    use_ED = FALSE){

  # Simplify SStype
  if(is.character(SStype))
    SStype <- nchar(SStype)
  if(!SStype %in% 1:3)
    stop("'SStype' must be one of 'I', 'II' or 'III'.")

  if(!missing(family))
    stop("hdanova() currently supports LM/LMM-style workflows only in this implementation.")
  if(any(grepl("|", formula, fixed = TRUE)))
    stop("hdanova() currently supports mixlm-style random effects r(), not lme4 '|' notation.")
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
  pca_active <- !(is.logical(pca.in) && !pca.in) && !identical(pca.in, 0)
  pls_active <- !(is.logical(pls.in) && !pls.in) && !identical(pls.in, 0)
  if(pca_active && pls_active)
    stop("Only one of 'pca.in' and 'pls.in' can be active at a time.")
  if(!is.null(REML) && !is.logical(REML))
    stop("'REML' must be NULL, TRUE or FALSE.")
  if(is.logical(REML) && !mixed_r){
    warning("'REML' was supplied but no random-effects r() terms were found in the formula; ignoring 'REML' and fitting a fixed-effects model.")
    REML <- NULL
  }
  ssq_method <- match.arg(REML_ssq_method, c("exact_refit", "wald", "ls"))
  ssq_method_kernel <- ssq_method
  use_ED_active <- isTRUE(use_ED) && mixed_r && is.logical(REML) && aug_is_numeric
  if(isTRUE(use_ED) && !use_ED_active)
    warning("'use_ED' currently only affects numeric 'aug_error' in REML/ML mixed-model fits; using standard degrees of freedom instead.")

  old_options <- options(contrasts = c(contrasts, "contr.poly"))
  on.exit(options(old_options), add = TRUE)

  ########################## Combined effects ##########################
  # Check for combined factors and remove symbols from formula.
  combined <- FALSE
  if(grepl("comb(", as.character(formula)[3], fixed = TRUE)){
    combined <- list()
    form_no_r <- mixlm::rparse(formula)
    tl <- attr(terms(form_no_r), "term.labels")
    j <- 1
    for(i in seq_along(tl)){
      if(grepl("comb(", tl[i], fixed = TRUE)){
        combined[[j]] <- attr(cparse(formula(paste0(".~", tl[i]))), "terms")[[1]]
        j <- j + 1
      }
    }
    formula <- cparse(formula)
  }

  ## Get the data matrices
  formula_mf <- formula
  if(mixed_r){
    rw <- mixlm::random.worker(formula, data, REML)
    formula_mf <- rw$formula
  }

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  if(missing(weights))
    mf$weights <- NULL
  if(missing(subset))
    mf$subset <- NULL
  if(missing(na.action))
    mf$na.action <- NULL
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
  if(pca.in != 0){ # Pre-decomposition, e.g., LiMM-PCA, PC-ANOVA
    # Automatic determination of dimension
    if(is.logical(pca.in) && pca.in)
      pca.in <- which.min(.PCAcv(Y))
    # Dimensions according to explained variance
    if(pca.in < 1){
      pca <- .pca(Y)
      pca.in <- min(which(cumsum(pca$explvar / 100) >= pca.in))
    }
    # Limit number of extracted components
    if(pca.in > p){
      warning(paste0("Reducing 'pca.in' from ", pca.in, " to the number of variables (", p, ")"))
      pca.in <- p
    }
    # PCA scores
    pca <- .pca(Y, ncomp = pca.in)
    Y <- pca$scores
  }

  ########################## PLS of response matrix ##########################
  if(pls.in != 0){
    formula_pls <- mixlm::rparse(formula)
    if(any(grepl("comb(", formula_pls, fixed = TRUE)))
      formula_pls <- cparse(formula_pls)

    tt_pls <- stats::delete.response(terms(formula_pls, data = mf_data))
    Yresp_pls <- model.matrix(tt_pls, data = mf_data)
    if("(Intercept)" %in% colnames(Yresp_pls))
      Yresp_pls <- Yresp_pls[, colnames(Yresp_pls) != "(Intercept)", drop = FALSE]
    if(ncol(Yresp_pls) == 0)
      stop("'pls.in' requires at least one non-constant predictor column after stripping r()/comb() and removing intercept.")

    max_pls <- min(pracma::Rank(Y), ncol(Y), nrow(Y) - 1, ncol(Yresp_pls))
    if(max_pls < 1)
      stop("'pls.in' has no estimable components for the current data.")

    dat_pls <- data.frame(X = I(Y), Y = I(Yresp_pls))
    if(is.logical(pls.in) && pls.in)
      pls.in <- max_pls
    if(!is.numeric(pls.in) || length(pls.in) != 1 || !is.finite(pls.in) || pls.in <= 0)
      stop("'pls.in' must be FALSE, TRUE, or a positive numeric value.")

    pls_fit <- NULL
    if(pls.in < 1){
      pls_fit <- pls::plsr(Y ~ X, data = dat_pls, ncomp = max_pls)
      pls.in <- which(cumsum(pls::explvar(pls_fit) / 100) >= pls.in)[1]
    }
    pls.in <- as.integer(round(pls.in))
    if(pls.in > max_pls){
      warning(paste0("Reducing 'pls.in' from ", pls.in, " to ", max_pls, "."))
      pls.in <- max_pls
    }
    if(is.null(pls_fit))
      pls_fit <- pls::plsr(Y ~ X, data = dat_pls, ncomp = pls.in)
    Y <- as.matrix(pls_fit$scores[, seq_len(pls.in), drop = FALSE])
  }

  mf_data[[1]] <- Y
  p_fit <- ncol(Y)
  fit.type <- if(mixed_r) "'lmm' (Linear Mixed Model)" else "'lm' (Linear Model)"

  weights_vec <- stats::model.weights(mf_data)
  if(!is.null(weights_vec)){
    weights_vec <- as.numeric(weights_vec)
    if(any(!is.finite(weights_vec)) || any(weights_vec < 0))
      stop("'weights' must be finite and non-negative.")
  }
  models <- NULL
  mom_anova <- NULL
  ########################## Decide model type ##########################
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
    if(is.null(colnames(Y)))
      colnames(Y) <- paste0("Y", seq_len(ncol(Y)))
    resp_name <- as.character(formula[[2]])
    models <- vector("list", p_fit)
    names(models) <- colnames(Y)
    for(i in seq_len(p_fit)){
      dat_i <- mf_data
      dat_i[[resp_name]] <- Y[, i, drop = FALSE]
      lm_args_i <- list(formula = formula,
                        data = dat_i,
                        unrestricted = unrestricted,
                        equal_baseline = equal_baseline,
                        contrasts = contrasts)
      if(!is.null(weights_vec))
        lm_args_i$weights <- weights_vec
      models[[i]] <- do.call(mixlm::lm, lm_args_i)
    }
    if(mixed_r)
      mom_anova <- mixlm::Anova(mod, type = c("I", "II", "III")[SStype])
  }

  ########################## Dry-run to find properties ##########################
  M <- model.matrix(mod)
  X_fixed <- M
  assign_fixed <- attr(M, "assign")
  assign <- assign_fixed
  modFra <- HDANOVA::extended.model.frame(model.frame(mod), mf_data)
  effs <- attr(terms(mod), "term.labels")
  effs_fixed <- effs

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

  maxDir <- numeric(n_terms)
  for(i in seq_len(n_terms))
    maxDir[i] <- min(sum(assign == i), p)
  effsAB <- effs

  ########################## ANOVA ##########################
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

  ########################## LS estimates ##########################
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

  ########################## Augmented LS/errors ##########################
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
  ssq_method_used <- "qr"
  ssq_source <- "qr"

  approved <- seq_len(n_terms)
  names(approved) <- effs
  approvedAB <- approved
  approvedComb <- as.list(approved)
  names(approvedComb) <- names(approvedAB)
  eff_combined <- setNames(rep(FALSE, n_terms), effs)

  if(is.logical(REML)){
    if(exists(".ML_variance_partition_all_types", mode = "function", inherits = TRUE)){
      vp <- .ML_variance_partition_all_types(models, type = SStype, method = ssq_method_kernel)
      ssq <- stats::setNames(rep(0, length(effs) + 1), c(effs, "Residuals"))
      common_terms <- intersect(names(ssq), vp$Term)
      if(length(common_terms) > 0)
        ssq[common_terms] <- vp$SSQ[match(common_terms, vp$Term)]
      ssq_method_used <- ssq_method
      ssq_source <- "all_types"
    } else if(exists(".ML_variance_partition_single", mode = "function", inherits = TRUE)){
      ssq <- stats::setNames(rep(0, length(effs) + 1), c(effs, "Residuals"))
      for(i in seq_along(models))
        ssq <- ssq + .ML_variance_partition_single(models[[i]], SStype)$Variance * N
      ssq_method_used <- "legacy_single"
      ssq_source <- "single"
    } else if("HDANOVA" %in% loadedNamespaces()){
      if(exists(".ML_variance_partition_all_types", envir = asNamespace("HDANOVA"), mode = "function", inherits = FALSE)){
        vp_fun_post <- get(".ML_variance_partition_all_types", envir = asNamespace("HDANOVA"))
        vp <- vp_fun_post(models, type = SStype, method = ssq_method_kernel)
        ssq <- stats::setNames(rep(0, length(effs) + 1), c(effs, "Residuals"))
        common_terms <- intersect(names(ssq), vp$Term)
        if(length(common_terms) > 0)
          ssq[common_terms] <- vp$SSQ[match(common_terms, vp$Term)]
        ssq_method_used <- ssq_method
        ssq_source <- "all_types"
      } else {
        vp_fun_post <- get(".ML_variance_partition_single", envir = asNamespace("HDANOVA"))
        ssq <- stats::setNames(rep(0, length(effs) + 1), c(effs, "Residuals"))
        for(i in seq_along(models))
          ssq <- ssq + vp_fun_post(models[[i]], SStype)$Variance * N
        ssq_method_used <- "legacy_single"
        ssq_source <- "single"
      }
    } else {
      warning("Variance-partition helper not found; using QR-based SSQ fallback for REML branch.")
      ssq_terms <- rowSums(SS_matrix)
      ssq <- c(stats::setNames(ssq_terms, effs), Residuals = sum(SSE_full))
      ssq_method_used <- "qr_fallback"
      ssq_source <- "qr"
    }
  } else {
    ssq_terms <- rowSums(SS_matrix)
    ssq <- c(stats::setNames(ssq_terms, effs), Residuals = sum(SSE_full))
    ssq_method_used <- if(isTRUE(respect_SStype)) "qr_sstype" else "qr_regression"
    ssq_source <- "qr"
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

  ########################## Combined effects ##########################
  remove <- integer(0)
  if(is.list(combined)){
    ## --- Validate comb() terms before processing ---
    for(i in seq_along(combined)){
      comb_terms <- combined[[i]]
      if(length(comb_terms) < 2)
        stop(paste0("'comb()' must contain at least two terms; got: comb(",
                    paste(comb_terms, collapse = " + "), ")."))
      for(dis in comb_terms){
        if(!(dis %in% effsAB)){
          # Check whether the term is a reordered interaction of a known effect
          if(grepl(":", dis, fixed = TRUE)){
            dis_sorted <- paste(sort(strsplit(dis, ":", fixed = TRUE)[[1]]), collapse = ":")
            effsAB_sorted <- vapply(effsAB, function(e)
              paste(sort(strsplit(e, ":", fixed = TRUE)[[1]]), collapse = ":"), character(1))
            match_perm <- effsAB[effsAB_sorted == dis_sorted]
            if(length(match_perm) > 0){
              stop(paste0("Term '", dis, "' not found in model effects. ",
                          "Interaction terms follow R's colon-order; ",
                          "did you mean '", match_perm[1], "'?"))
            }
          }
          stop(paste0("Term '", dis, "' specified in comb() was not found among model effects: ",
                      paste(effsAB, collapse = ", "), "."))
        }
      }
    }
    ## --- end validation ---
    for(i in seq_along(combined)){
      combName <- paste(combined[[i]], collapse = "+")
      effs <- c(effs, combName)
      maxDir <- c(maxDir, 0L)
      approved <- c(approved, length(effs))
      approvedAB <- c(approvedAB, length(effs))
      approvedComb <- c(approvedComb, list(integer(0)))
      names(approved)[length(approved)] <- combName
      names(approvedAB)[length(approvedAB)] <- combName

      LS[[approved[length(approved)]]] <- 0 * LS[[effs[[approved[1]]]]]
      names(LS)[approved[length(approved)]] <- combName
      LS_regression_out[[approved[length(approved)]]] <- 0 * LS_regression_out[[effs[[approved[1]]]]]
      names(LS_regression_out)[approved[length(approved)]] <- combName
      LS_SStype_out[[approved[length(approved)]]] <- 0 * LS_SStype_out[[effs[[approved[1]]]]]
      names(LS_SStype_out)[approved[length(approved)]] <- combName
      LS_aug[[approved[length(approved)]]] <- 0 * LS_aug[[effs[[approved[1]]]]]
      names(LS_aug)[approved[length(approved)]] <- combName
      aug_regression$LS_aug[[approved[length(approved)]]] <- 0 * aug_regression$LS_aug[[effs[[approved[1]]]]]
      names(aug_regression$LS_aug)[approved[length(approved)]] <- combName
      aug_SStype$LS_aug[[approved[length(approved)]]] <- 0 * aug_SStype$LS_aug[[effs[[approved[1]]]]]
      names(aug_SStype$LS_aug)[approved[length(approved)]] <- combName
      eff_combined[approved[length(approved)]] <- TRUE
      names(eff_combined)[approved[length(approved)]] <- combName
      error[[approved[length(approved)]]] <- 0 * residuals
      names(error)[approved[length(approved)]] <- combName
      aug_regression$error[[approved[length(approved)]]] <- 0 * residuals
      names(aug_regression$error)[approved[length(approved)]] <- combName
      aug_SStype$error[[approved[length(approved)]]] <- 0 * residuals
      names(aug_SStype$error)[approved[length(approved)]] <- combName
      ssq[approved[length(approved)]] <- 0
      names(ssq)[approved[length(approved)]] <- combName
      maxDir[approved[length(approved)]] <- 0

      for(dis in combined[[i]]){
        if(any(!(dis %in% names(approvedAB))))
          stop("Cannot combine a continuous effect with a categorical factor.")
        remove <- c(remove, which(effsAB == dis))
        idx_dis <- approvedAB[names(approvedAB) == dis]
        LS[[approved[length(approved)]]] <- LS[[approved[length(approved)]]] + LS[[effs[idx_dis]]]
        LS_regression_out[[approved[length(approved)]]] <- LS_regression_out[[approved[length(approved)]]] +
          LS_regression_out[[effs[idx_dis]]]
        LS_SStype_out[[approved[length(approved)]]] <- LS_SStype_out[[approved[length(approved)]]] +
          LS_SStype_out[[effs[idx_dis]]]
        LS_aug[[approved[length(approved)]]] <- LS_aug[[approved[length(approved)]]] + LS_aug[[effs[idx_dis]]]
        aug_regression$LS_aug[[approved[length(approved)]]] <- aug_regression$LS_aug[[approved[length(approved)]]] +
          aug_regression$LS_aug[[effs[idx_dis]]]
        aug_SStype$LS_aug[[approved[length(approved)]]] <- aug_SStype$LS_aug[[approved[length(approved)]]] +
          aug_SStype$LS_aug[[effs[idx_dis]]]
        ssq[approved[length(approved)]] <- ssq[approved[length(approved)]] + ssq[effs[idx_dis]]
        maxDir[approved[length(approved)]] <- max(maxDir[approved[length(approved)]], maxDir[idx_dis])
        approvedComb[[length(approvedComb)]] <- c(approvedComb[[length(approvedComb)]], idx_dis)
        approved <- approved[names(approvedAB) != dis]
        approvedAB <- approvedAB[names(approvedAB) != dis]
      }

      error[[approved[length(approved)]]] <- LS[[approved[length(approved)]]] + residuals
      aug_regression$error[[approved[length(approved)]]] <- LS_regression_out[[approved[length(approved)]]] + residuals
      aug_SStype$error[[approved[length(approved)]]] <- LS_SStype_out[[approved[length(approved)]]] + residuals
      LS_aug[[approved[length(approved)]]] <- error[[approved[length(approved)]]]
      aug_regression$LS_aug[[approved[length(approved)]]] <- aug_regression$error[[approved[length(approved)]]]
      aug_SStype$LS_aug[[approved[length(approved)]]] <- aug_SStype$error[[approved[length(approved)]]]
      names(approvedComb)[length(approvedComb)] <- combName
    }

    if(names(ssq)[length(ssq)] == "Residuals")
      ssq <- c(ssq[setdiff(seq_len(length(ssq) - 1), remove)], ssq[length(ssq)])
    else
      ssq <- c(ssq[setdiff(seq_along(ssq), remove)], Residuals = sum(SSE_full))
    names(ssq)[length(ssq)] <- "Residuals"
    eff_combined <- eff_combined[approved]
  }

  model <- model.frame(mod)
  model[[1]] <- Yorig

  ########################## Return ##########################
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
                          assign_fixed = assign_fixed,
                          X_fixed = X_fixed,
                          effs_fixed = effs_fixed,
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
                          pls.in = pls.in,
                          Y_fit = Y,
                          weights = weights_vec,
                          respect_SStype = respect_SStype,
                          effect_source = effect_source,
                          ssq_method = ssq_method_used,
                          ssq_source = ssq_source,
                          version = 2))
  if(!is.null(EDall))
    obj$ED <- EDall
  if(pca.in != 0)
    obj$Ypca <- list(pca = pca, ncomp = pca.in)
  if(pls.in != 0)
    obj$Ypls <- list(pls = pls_fit, ncomp = pls.in)
  class(obj) <- c("hdanova", "list")
  obj
}

dummyvar <- function(x){
  dv <- model.matrix(~x-1, data.frame(x=x, check.names=FALSE))
  colnames(dv) <- levels(x)
  dv
}

SS <- function(x) sum(x^2)

# Not used
formula_comb <- function(formula, data, REML = NULL){
  formula <- formula(formula)
  terms <- terms(formula)
  effsr <- attr(terms,"term.labels")
#  effs  <- attr(terms(cparse(formula)),"term.labels")
  effs  <- attr(cparse(formula),"terms")[[1]]
  if(length(effs)==0){
    return( list(0) )
  }

  has.intercept <- attr(terms,"intercept")==1
  combs <- sort(unique(c(grep("[:]comb[(]",effsr),   # Match combination in interaction
                         grep("^comb[(]",  effsr),   # Match combination in the beginning
                         grep("[(]comb[(]",effsr))))  # Match combination inside function

  eff.splits <- list()
  for(i in 1:length(effs)){ # Split effect to look for hidden combination interactions
    eff.splits[[i]] <- fparse(formula(paste("1~", effs[i],sep="")))
  }
  eff.lengths <- lapply(eff.splits,length)
  main.effs   <- effs[eff.lengths==1]
  main.combs  <- main.effs[main.effs%in%effs[combs]]
  main.combs.only.inter <- character(0)
  for(i in combs){
    main.combs.only.inter <- c(main.combs.only.inter, setdiff(eff.splits[[i]],main.effs)) # Combination main effects only present in interactions
  }
  inter.combs <- which(unlist(lapply(eff.splits,function(i) any(main.combs%in%i))))
  # Check if any interactions containing combination effects are not labeled as combination
  if(any(is.na(match(inter.combs,combs)))){
    extra.randoms <- inter.combs[which(is.na(match(inter.combs,combs)))]
    warning(paste(paste(effs[extra.randoms],sep="",collapse=", "), " included as combination interaction",ifelse(length(extra.randoms)==1,"","s"),sep=""))
    combs <- cbind(combs,extra.randoms)
    effs  <- effs[!(extra.randoms%in%effs)]
  }
  if(length(combs)==0){
    return( list(0) )
  } else {
    if(is.logical(REML)){
      remleffs     <- c(effs[setdiff(1:length(effs),combs)],paste("(1|",effs[combs],")",sep=""))
      reml.formula <- formula(paste(formula[[2]],"~",paste(remleffs,collapse="+"),ifelse(has.intercept,"","-1"),sep=""))

      return( list(cformula = formula,
                   formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")),
                   combination = effs[combs],
                   main.combs.only.inter = main.combs.only.inter,
                   fixed = effs[setdiff(1:length(effsr),combs)],
                   all = effs,
                   allc = effsr,
                   has.intercept = has.intercept,
                   remleffs = remleffs,
                   reml.formula = reml.formula))
    } else {
      return( list(cformula = formula,
                   formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")),
                   combination = effs[combs],
                   main.combs.only.inter = main.combs.only.inter,
                   fixed = effs[setdiff(1:length(effsr),combs)],
                   all = effs,
                   allc = effsr,
                   has.intercept = has.intercept))
    }
  }
}

# Remove comb() from formula (possibly convert to lmer formula)
cparse <- function (f, REML = NULL) {
  if(!inherits(f,'formula'))
    # if (class(f) != "formula")
    stop("'f' must be a formula")
  right <- attr(terms(f),"term.labels") # Let R split short-hand notations first
  if( length(right) == 0){
    return(f)
  }
  n <- length(right)
  result      <- character(n)
  result.REML <- character(n)
  terms       <- list()

  # Main recursive loop extracting effects without comb()
  for(i in 1:n){
    parsecall <- function(x) {
      if(!inherits(x,'call'))
        # if (class(x) != "call")
        stop("'x' must be a call")
      if (length(x[[1]]) == 1) {
        return(x)
      }
      if (length(x[[1]]) == 2 && x[[1]][[1]] == "comb") {
        return(x[[1]][2])
      }
      for (i in 2:length(x[[1]])) x[[1]][i] <- parsecall(x[[1]][i])
      return(x)
    }
    result[[i]] <- as.character(parsecall(formula(paste("~",paste(right[i],sep="+")))[2]))
    terms[[i]]  <- trimws(strsplit(result[[i]], "\\+")[[1]])
  }
  f[3] <- formula(paste("~", paste(result, sep="", collapse="+")))[2]

  # Recursive loop adding (1 | x) notation for REML estimation
  if(is.logical(REML)){
    for(i in 1:n){
      parsecall <- function(x) {
        if(!inherits(x,'call'))
          # if (class(x) != "call")
          stop("'x' must be a call")
        if (length(x[[1]]) == 1) {
          return(FALSE)
        }
        if (length(x[[1]]) == 2 && x[[1]][[1]] == "comb") {
          return(TRUE)
        } else {
          tmp <- logical(length(x[[1]]))
          for (j in 2:length(x[[1]])) tmp[j-1] <- parsecall(x[[1]][j])
          return(any(tmp))
        }
      }
      ran <- parsecall(formula(paste("~",paste(right[i],sep="+")))[2])
      result.REML[i] <- ifelse(ran,as.character(formula(paste("~(1 | ",result[[i]],")",sep=""))[2]), result[[i]])
      terms[[i]]  <- trimws(strsplit(result.REML[[i]], "\\+")[[1]])
    }
    f[3] <- formula(paste("~", paste(result.REML, sep="", collapse="+")))[2]
  }
  attr(f, "terms") <- terms
  f
}

# Extract variable names from formula
fparse <- function(f) {
  if(!inherits(f,'formula')) stop("'f' must be a formula")
  # if (class(f) != "formula") stop("'f' must be a formula")
  right <- f[3]
  parsecall <- function(x) {
    if(!inherits(x,'call')) stop("'x' must be a call")
    # if (class(x) != "call") stop("'x' must be a call")
    if (length(x[[1]]) == 1) {
      if(is.numeric(x[[1]])) {
        return()
      } else {
        return(as.character(x[[1]]))
      }
    }
    res <- list()
    for (i in 2:length(x[[1]]))
      res[[i-1]] <- parsecall(x[[1]][i])
    return(unlist(res))
  }
  unique(parsecall(right))
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

  ED <- vector("list", ncol(Y))
  N <- nrow(Y)
  q <- ncol(Z)

  for(i in seq_len(ncol(Y))){
    wrong.ED <- FALSE
    G <- diag(q) * 0
    r0 <- 0
    for(j in seq_along(randEffList)){
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
