#' Permutation for SS-Type-Aligned HDANOVA
#'
#' @description
#' Permutation testing for \code{hdanova()} objects using SS-type-aligned QR
#' effect matrices. This function is intended as an opt-in alternative to the
#' legacy regression-based \code{permutation()} workflow.
#'
#' @details
#' The permutation statistics are computed from SS-type-aligned QR reduced/full
#' model contrasts rather than the legacy regression-based LS matrices. Fixed
#' models and mixed MoM models are supported. Approximate permutation uses a
#' relaxed global shuffle of observations; exact permutation uses permissible
#' block-restricted shuffles. For REML/ML mixed models with
#' \code{respect_SStype = TRUE}, a warning is issued to highlight that
#' permutation statistics and fitted-model REML SSQ decompositions are based on
#' different computational definitions.
#'
#' @param object A \code{hdanova} object.
#' @param permute Number of permutations to perform (default = 1000).
#' @param perm.type Type of permutation to perform, either \code{"approximate"}
#' or \code{"exact"} (default = \code{"approximate"}).
#' @param respect_SStype Logical or \code{NULL}. If \code{NULL} (default),
#' use the \code{hdanova} object setting (\code{object$more$respect_SStype}).
#' If \code{FALSE}, follow the legacy regression-based permutation logic. If
#' \code{TRUE}, use SS-type-aligned QR contrasts. For REML/ML mixed models,
#' this may yield permutation statistics that differ from fitted-model REML
#' SSQ decompositions selected by \code{REML_ssq_method}.
#' @param unique.digits Number of digits used when rounding permutation SSQ values before
#' checking uniqueness (default = 12). Set to \code{NULL} to disable this warning.
#' @param unique.frac Minimum fraction of unique rounded SSQ values required to avoid
#' warning (default = 0.95). Set to \code{NULL} to disable this warning.
#' @param exhaustive.warn Logical; if \code{TRUE} (default), warn when exact permutation
#' uses exhaustive enumeration with fewer permutations than requested.
#'
#' @returns An updated \code{hdanova} object with \code{permute} results.
#'
#' @export
permutation <- function(object,
                         permute = 1000,
                         perm.type = c("approximate", "exact"),
                         respect_SStype = NULL,
                         unique.digits = 12,
                         unique.frac = 0.95,
                         exhaustive.warn = TRUE){
  perm.type <- match.arg(perm.type)
  if(!inherits(object, "hdanova"))
    stop("'object' must inherit from 'hdanova'.")
  if(is.null(respect_SStype)){
    respect_SStype <- isTRUE(object$more$respect_SStype)
  } else {
    if(!is.logical(respect_SStype) || length(respect_SStype) != 1 || is.na(respect_SStype))
      stop("'respect_SStype' must be TRUE, FALSE, or NULL.")
  }

  if(isTRUE(is.logical(object$more$REML)) && isTRUE(respect_SStype) &&
     !is.null(object$more$ssq_method) &&
      object$more$ssq_method %in% c("exact_refit", "wald", "ls")){
    warning("REML/ML permutation with 'respect_SStype = TRUE' uses SS-type-aligned QR permutation statistics. These may differ from fitted-model REML SSQ decomposition methods ('exact_refit', 'wald', 'ls').")
  }

  # Keep parity with hdanova default behavior: regression path unless explicitly aligned.
  if(!respect_SStype){
    object_reg <- .permutation_regression_impl(
      object = object,
      permute = permute,
      perm.type = perm.type,
      unique.digits = unique.digits,
      unique.frac = unique.frac,
      exhaustive.warn = exhaustive.warn
    )
    return(object_reg)
  }

  check_unique_warning <- !is.null(unique.digits) && !is.null(unique.frac)
  if(check_unique_warning){
    if(!is.numeric(unique.digits) || length(unique.digits) != 1 || is.na(unique.digits))
      stop("'unique.digits' must be a single numeric value or NULL.")
    if(!is.numeric(unique.frac) || length(unique.frac) != 1 || is.na(unique.frac) ||
       unique.frac < 0 || unique.frac > 1)
      stop("'unique.frac' must be a single numeric value in [0, 1] or NULL.")
    unique.digits <- as.integer(unique.digits)
  }
  if(!is.logical(exhaustive.warn) || length(exhaustive.warn) != 1 || is.na(exhaustive.warn))
    stop("'exhaustive.warn' must be TRUE or FALSE.")

  approved <- object$more$approved
  eff_names <- object$more$effs[approved]
  Y_model <- object$more$Y_fit
  if(is.null(Y_model))
    stop("Object does not contain the transformed response matrix required by permutation().")
  sqrt_weights <- if(is.null(object$more$weights)) rep(1, object$more$N) else sqrt(object$more$weights)

  ssqa <- pvals <- numeric(length(approved))
  ssqaperm <- vector("list", length(approved))
  permutation_orders <- vector("list", length(approved))
  exchangeable_units <- list()
  perm_warnings <- character(0)
  names(ssqa) <- names(pvals) <- names(ssqaperm) <- names(permutation_orders) <- eff_names

  caches <- .hda_build_permutation_caches(object, sqrt_weights)

  if(!(is.logical(permute) && !permute)){
    if(is.logical(permute))
      permute <- 1000

    for(i in seq_along(approved)){
      eff_name <- eff_names[i]
      cache <- caches[[eff_name]]
      if(is.null(cache)){
        # Random effect term: not testable via fixed-effects permutation; skip.
        ssqa[eff_name] <- NA_real_
        pvals[eff_name] <- NA_real_
        ssqaperm[[eff_name]] <- numeric(0)
        permutation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
        next
      }
      ssqa[eff_name] <- .hda_effect_ss_from_cache(cache, Y_model, sqrt_weights)

      effect_vars <- .effect_variable_names(object, eff_name)
      block_info <- .build_permutation_blocks(object, eff_name, effect_vars)
      perm_engine <- .hda_prepare_effect_permutation_engine(
        object = object,
        cache = cache,
        effect_vars = effect_vars,
        Y = Y_model,
        sqrt_weights = sqrt_weights
      )
      exchangeable_units[[eff_name]] <- list(
        total_blocks = length(block_info$groups),
        exchangeable_blocks = block_info$multi_units,
        block_sizes = lengths(block_info$groups)
      )

      if(perm.type == "exact" && block_info$multi_units == 0){
        perm_warnings <- c(perm_warnings,
                           paste0("Effect '", eff_name, "' has no exchangeable units; permutation skipped."))
        ssqaperm[[eff_name]] <- numeric(0)
        permutation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
        pvals[eff_name] <- NA_real_
        next
      }

      if(perm.type == "exact"){
        n_possible <- block_info$n_permutations
        if(isTRUE(exhaustive.warn) && is.finite(n_possible) && n_possible < permute){
          perm_warnings <- c(perm_warnings,
                             paste0("Effect '", eff_name, "' uses exhaustive permutation with ",
                                    .format_number(n_possible),
                                    " unique orders (requested ", .format_number(permute), ")."))
        }
        if(is.finite(n_possible) && n_possible <= permute){
          perms <- numeric(n_possible)
          order_matrix <- matrix(0L, nrow = n_possible, ncol = object$more$N)
          counter <- 0L
          group_perms <- lapply(block_info$groups, .all_permutations)
          .enumerate_permutation_indices(block_info$groups, group_perms, function(order){
            counter <<- counter + 1L
            order_matrix[counter, ] <<- order
            perms[counter] <<- .hda_effect_ss_from_permuted_effect(
              engine = perm_engine,
              perm_index = order,
              Y = Y_model,
              sqrt_weights = sqrt_weights
            )
          })
          permutation_orders[[eff_name]] <- order_matrix
        } else {
          perms <- numeric(permute)
          pb <- progress::progress_bar$new(total = permute,
                                           format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)"))
          for(perm in seq_len(permute)){
            idx <- .permute_within_blocks(block_info$groups)
            perms[perm] <- .hda_effect_ss_from_permuted_effect(
              engine = perm_engine,
              perm_index = idx,
              Y = Y_model,
              sqrt_weights = sqrt_weights
            )
            pb$tick()
          }
          permutation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
        }
      } else {
        perms <- numeric(permute)
        pb <- progress::progress_bar$new(total = permute,
                 format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)"))
        for(perm in seq_len(permute)){
          # Relaxed approximate permutation: global shuffle of the effect labels.
          idx <- sample(object$more$N)
          perms[perm] <- .hda_effect_ss_from_permuted_effect(
            engine = perm_engine,
            perm_index = idx,
            Y = Y_model,
            sqrt_weights = sqrt_weights
          )
          pb$tick()
        }
        permutation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
      }

      ssqaperm[[eff_name]] <- perms
      if(check_unique_warning && length(perms) > 0){
        unique_ssq <- length(unique(round(perms, unique.digits)))
        if(unique_ssq < length(perms) * unique.frac){
          perm_warnings <- c(perm_warnings,
                             paste0("Effect '", eff_name, "' produced only ",
                                    .format_number(unique_ssq),
                                    " unique SSQ values (rounded to ", unique.digits, " decimals) out of ",
                                    .format_number(length(perms)),
                                    " evaluated permutations."))
        }
      }
      pvals[eff_name] <- sum(perms > ssqa[eff_name]) / length(perms)
    }
  }

  if(length(perm_warnings) > 0)
    warning(paste(unique(perm_warnings), collapse = "\n"))

  object$permute <- list(ssqa = ssqa,
                          ssqaperm = ssqaperm,
                          pvalues = pvals,
                          permutations = permute,
                          exchangeable = exchangeable_units,
                          orders = permutation_orders,
                          method = "permutation",
                          perm.type = perm.type,
                          effect_source = "SStype",
                          ssq_method = object$more$ssq_method)
  object
}


.hda_build_permutation_caches <- function(object, sqrt_weights){
  # For REML models, object$X includes random-effect Z columns which must NOT be
  # used in permutation caches — including Z in the nuisance model absorbs the
  # between-subject variation (including between-subject fixed effects like
  # Treatment), leaving nothing to permute. Use only fixed-effect columns.
  if(isTRUE(is.logical(object$more$REML)) && !is.null(object$more$X_fixed)){
    X_perm   <- object$more$X_fixed
    assign_p <- object$more$assign_fixed
    effs_p   <- object$more$effs_fixed
  } else {
    X_perm   <- object$X
    assign_p <- object$more$assign
    effs_p   <- object$more$effs
  }

  effs <- object$more$effs
  caches <- vector("list", length(object$more$approved))
  names(caches) <- effs[object$more$approved]
  for(i in seq_along(object$more$approved)){
    eff_name <- effs[object$more$approved[i]]
    # Map approved index to the position within fixed-only effs for column lookup
    eff_pos_fixed <- match(eff_name, effs_p)
    if(is.na(eff_pos_fixed)){
      # Random effect term — skip building a permutation cache for it
      next
    }
    cols <- .hda_effect_column_sets(assign_p, effs_p, eff_pos_fixed, object$SStype)
    caches[[eff_name]] <- list(
      curr = .hda_qr_cache(X_perm[, cols$curr, drop = FALSE], sqrt_weights),
      prev = .hda_qr_cache(X_perm[, cols$prev, drop = FALSE], sqrt_weights)
    )
  }
  caches
}


.hda_effect_column_sets <- function(assign, term_labels, term_index, SStype){
  if(SStype == 1){
    curr <- which(assign <= term_index)
    prev <- which(assign < term_index)
  } else if(SStype == 2){
    is_higher <- grepl(paste0("^", term_labels[term_index], ":"), term_labels) |
      grepl(paste0(":", term_labels[term_index], "$"), term_labels)
    curr <- which(!(assign %in% which(is_higher)))
    prev <- which(!(assign %in% c(term_index, which(is_higher))))
  } else {
    curr <- seq_along(assign)
    prev <- which(assign != term_index)
  }
  list(curr = curr, prev = prev)
}


.hda_qr_cache <- function(X, sqrt_weights){
  list(X = X,
       QR = if(ncol(X) > 0) qr(X * sqrt_weights) else NULL)
}


.hda_weighted_fitted_from_cache <- function(cache, Y, sqrt_weights){
  fitted <- matrix(0, nrow = nrow(Y), ncol = ncol(Y), dimnames = dimnames(Y))
  if(is.null(cache$QR) || ncol(cache$X) == 0)
    return(fitted)
  coefs <- qr.coef(cache$QR, Y * sqrt_weights)
  estimable <- !is.na(rowSums(coefs))
  if(any(estimable))
    fitted <- cache$X[, estimable, drop = FALSE] %*% coefs[estimable, , drop = FALSE]
  fitted
}


.hda_effect_ss_from_cache <- function(cache, Y, sqrt_weights){
  Yhat_curr <- .hda_weighted_fitted_from_cache(cache$curr, Y, sqrt_weights)
  Yhat_prev <- .hda_weighted_fitted_from_cache(cache$prev, Y, sqrt_weights)
  sum((sqrt_weights * (Yhat_curr - Yhat_prev))^2)
}


.hda_prepare_effect_permutation_engine <- function(object, cache, effect_vars, Y, sqrt_weights){
  model_df <- object$model.frame
  if(is.null(model_df))
    model_df <- model.frame(object$models[[1]])
  if(ncol(model_df) > 0)
    model_df <- model_df[, -1, drop = FALSE]

  x_target <- colnames(object$X)
  curr_names <- colnames(cache$curr$X)
  prev_names <- colnames(cache$prev$X)
  Yhat_prev <- .hda_weighted_fitted_from_cache(cache$prev, Y, sqrt_weights)
  resid_prev <- Y - Yhat_prev

  list(model_df = model_df,
       effect_vars = intersect(effect_vars, names(model_df)),
       x_target = x_target,
       curr_names = curr_names,
       prev_names = prev_names,
       cache = cache,
       Yhat_prev = Yhat_prev,
       resid_prev = resid_prev)
}


.hda_align_design_columns <- function(X, target_names){
  if(length(target_names) == 0)
    return(matrix(0, nrow = nrow(X), ncol = 0))

  out <- matrix(0, nrow = nrow(X), ncol = length(target_names))
  colnames(out) <- target_names
  nx <- .hda_normalize_design_names(colnames(X))
  nt <- .hda_normalize_design_names(target_names)
  map <- match(nt, nx)
  ok <- which(!is.na(map))
  if(length(ok) > 0)
    out[, ok] <- X[, map[ok], drop = FALSE]
  out
}


.hda_normalize_design_names <- function(nms){
  gsub("[()]", "", nms)
}


.hda_effect_ss_from_permuted_effect <- function(engine, perm_index, Y, sqrt_weights){
  # Freedman-Lane style permutation for SS-type tests:
  # keep nuisance fit fixed and permute reduced-model residuals.
  Y_perm <- engine$Yhat_prev + engine$resid_prev[perm_index, , drop = FALSE]
  .hda_effect_ss_from_cache(engine$cache, Y_perm, sqrt_weights)
}


.permutation_regression_impl <- function(object,
                                          permute,
                                          perm.type,
                                          unique.digits,
                                          unique.frac,
                                          exhaustive.warn){
  check_unique_warning <- !is.null(unique.digits) && !is.null(unique.frac)
  if(check_unique_warning){
    if(!is.numeric(unique.digits) || length(unique.digits) != 1 || is.na(unique.digits))
      stop("'unique.digits' must be a single numeric value or NULL.")
    if(!is.numeric(unique.frac) || length(unique.frac) != 1 || is.na(unique.frac) ||
       unique.frac < 0 || unique.frac > 1)
      stop("'unique.frac' must be a single numeric value in [0, 1] or NULL.")
    unique.digits <- as.integer(unique.digits)
  }
  if(!is.logical(exhaustive.warn) || length(exhaustive.warn) != 1 || is.na(exhaustive.warn))
    stop("'exhaustive.warn' must be TRUE or FALSE.")

  ssqa <- pvals <- numeric(0)
  ssqaperm <- list()
  permutation_orders <- list()
  perm_warnings <- character(0)
  exchangeable_units <- list()

  if(!(is.logical(permute) && !permute)){
    if(is.logical(permute)){
      permute <- 1000
      if(interactive())
        message("Defaulting to 1000 permutations\n")
    }

    approved <- object$more$approved
    eff_names <- object$more$effs[approved]
    ssqa <- pvals <- numeric(length(approved))
    ssqaperm <- vector("list", length(approved))
    permutation_orders <- vector("list", length(approved))
    names(ssqa) <- names(pvals) <- names(ssqaperm) <- names(permutation_orders) <- eff_names
    
    # Check if C++ permutation kernel is available
    cpp_available <- exists(".permutation_ssq_kernel_cpp", mode = "function")
    cpp_failed <- FALSE

    if(perm.type == "approximate"){
      for(i in seq_along(approved)){
        a <- object$more$approvedAB[i]
        eff_name <- object$more$effs[a]
        D <- object$X[, object$more$assign %in% object$more$approvedComb[[names(a)]], drop = FALSE]
        DD <- D %*% pracma::pinv(D)
        DR <- object$more$LS_aug[[eff_name]]
        ssqa[eff_name] <- norm(DD %*% DR, "F")^2

        perms <- NULL
        if(cpp_available && !cpp_failed){
          cpp_res <- tryCatch(
            list(values = .run_cpp_permutation_with_progress(
              DD = DD,
              DR = DR,
              n_perm = permute,
              format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)")
            )),
            error = function(e)
              list(values = NULL, error = e)
          )
          perms <- cpp_res$values
          if(!is.null(cpp_res$error)){
            cpp_failed <- TRUE
            warning(paste0("C++ permutation kernel failed for effect '", eff_name,
                          "'; falling back to R implementation. ", cpp_res$error$message))
          }
        }
        
        if(is.null(perms)){
          perms <- numeric(permute)
          pb <- progress::progress_bar$new(total = permute,
                                           format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)"))
          for(perm in seq_len(permute)){
            perms[perm] <- norm(DD %*% DR[sample(object$more$N), ], "F")^2
            pb$tick()
          }
        }
        
        ssqaperm[[eff_name]] <- perms
        permutation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)

        if(check_unique_warning && length(perms) > 0){
          unique_ssq <- length(unique(round(perms, unique.digits)))
          if(unique_ssq < length(perms) * unique.frac){
            perm_warnings <- c(perm_warnings,
                               paste0("Effect '", eff_name, "' produced only ",
                                      .format_number(unique_ssq),
                                      " unique SSQ values (rounded to ", unique.digits, " decimals) out of ",
                                      .format_number(length(perms)),
                                      " evaluated permutations."))
          }
        }
        pvals[eff_name] <- sum(perms > ssqa[eff_name]) / permute
      }
    } else {
      for(i in seq_along(approved)){
        a <- object$more$approvedAB[i]
        eff_name <- object$more$effs[a]
        D <- object$X[, object$more$assign %in% object$more$approvedComb[[names(a)]], drop = FALSE]
        DD <- D %*% pracma::pinv(D)
        DR <- object$more$LS_aug[[eff_name]]
        ssqa[eff_name] <- norm(DD %*% DR, "F")^2

        effect_vars <- .effect_variable_names(object, eff_name)
        block_info <- .build_permutation_blocks(object, eff_name, effect_vars)
        exchangeable_units[[eff_name]] <- list(
          total_blocks = length(block_info$groups),
          exchangeable_blocks = block_info$multi_units,
          block_sizes = lengths(block_info$groups)
        )

        if(block_info$multi_units == 0){
          perm_warnings <- c(perm_warnings,
                             paste0("Effect '", eff_name, "' has no exchangeable units; exact permutation skipped."))
          ssqaperm[[eff_name]] <- numeric(0)
          permutation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
          pvals[eff_name] <- NA_real_
          next
        }

        n_possible <- block_info$n_permutations
        exhaustive_used <- is.finite(n_possible) && n_possible <= permute
        if(isTRUE(exhaustive.warn) && is.finite(n_possible) && n_possible < permute){
          perm_warnings <- c(perm_warnings,
                             paste0("Effect '", eff_name, "' uses exhaustive permutation with ",
                                    .format_number(n_possible),
                                    " unique orders (requested ", .format_number(permute), ")."))
        }

        if(is.finite(n_possible) && n_possible <= permute){
          perms <- numeric(n_possible)
          order_matrix <- matrix(0L, nrow = n_possible, ncol = object$more$N)
          counter <- 0L
          group_perms <- lapply(block_info$groups, .all_permutations)
          .enumerate_permutation_indices(block_info$groups, group_perms, function(order){
            counter <<- counter + 1L
            order_matrix[counter, ] <<- order
            perms[counter] <<- norm(DD %*% DR[order, ], "F")^2
          })
          permutation_orders[[eff_name]] <- order_matrix
        } else {
          perms <- numeric(permute)
          pb <- progress::progress_bar$new(total = permute,
                                           format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)"))
          for(perm in seq_len(permute)){
            idx <- .permute_within_blocks(block_info$groups)
            perms[perm] <- norm(DD %*% DR[idx, ], "F")^2
            pb$tick()
          }
          permutation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
        }

        ssqaperm[[eff_name]] <- perms
        if(check_unique_warning && length(perms) > 0){
          unique_ssq <- length(unique(round(perms, unique.digits)))
          if(unique_ssq < length(perms) * unique.frac){
            if(exhaustive_used){
              perm_warnings <- c(perm_warnings,
                                 paste0("Effect '", eff_name, "' produced only ",
                                        .format_number(unique_ssq),
                                        " unique SSQ values (rounded to ", unique.digits,
                                        " decimals) in exhaustive permutation."))
            } else {
              perm_warnings <- c(perm_warnings,
                                 paste0("Effect '", eff_name, "' produced only ",
                                        .format_number(unique_ssq),
                                        " unique SSQ values (rounded to ", unique.digits, " decimals) out of ",
                                        .format_number(length(perms)),
                                        " evaluated permutations."))
            }
          }
        }
        pvals[eff_name] <- sum(perms > ssqa[eff_name]) / length(perms)
      }
    }
  }

  if(length(perm_warnings) > 0)
    warning(paste(unique(perm_warnings), collapse = "\n"))

  object$permute <- list(ssqa = ssqa,
                          ssqaperm = ssqaperm,
                          pvalues = pvals,
                          permutations = permute,
                          exchangeable = exchangeable_units,
                          orders = permutation_orders,
                          method = "permutation",
                          perm.type = perm.type,
                          effect_source = "regression",
                          ssq_method = object$more$ssq_method)
  object
}


.effect_variable_names <- function(object, eff_name){
  # Return all component variable names for this effect, whether it is a main
  # effect or an interaction. For interactions (e.g. "A:B"), the component
  # variables ("A", "B") are used by .build_permutation_blocks to *exclude*
  # factor levels from the strata, leaving only true blocking factors (e.g.
  # random-effect grouping variables like Subject) to define exchange units.
  # This ensures:
  #   - Fixed model interaction: no blocking strata → global permutation.
  #   - Mixed model interaction: Subject blocks → within-subject permutation.
  raw_terms <- strsplit(eff_name, "+", fixed = TRUE)[[1]]
  raw_terms <- trimws(raw_terms)
  perm_vars <- character(0)
  for(term in raw_terms){
    parts <- .parse_effect_terms(term)
    if(length(parts) == 0)
      next
    perm_vars <- c(perm_vars, parts)
  }
  unique(perm_vars)
}


.parse_effect_terms <- function(name){
  if(is.null(name))
    return(character(0))
  clean <- gsub("\\s+", "", name)
  if(clean == "" || clean == "(Intercept)")
    return(character(0))
  terms <- strsplit(clean, "[:+\\|]", perl = TRUE)[[1]]
  terms[terms != ""]
}


.build_permutation_blocks <- function(object, eff_name, effect_vars){
  mf <- object$model.frame
  if(ncol(mf) > 0)
    mf <- mf[setdiff(names(mf), names(mf)[1])]
  candidate <- intersect(names(mf), object$more$effs)
  candidate <- setdiff(candidate, eff_name)
  candidate <- candidate[!candidate %in% "(Intercept)"]
  eligible <- vapply(candidate, function(name){
    vars <- .parse_effect_terms(name)
    length(intersect(vars, effect_vars)) == 0
  }, FUN.VALUE = TRUE)
  strata <- candidate[eligible]
  if(length(strata) == 0){
    keys <- rep("1", object$more$N)
  } else {
    data <- mf[strata]
    data[] <- lapply(data, function(col){
      col_chr <- as.character(col)
      col_chr[is.na(col_chr)] <- "NA"
      col_chr
    })
    keys <- do.call(paste, c(data, sep = "\r"))
  }
  groups <- split(seq_len(object$more$N), keys)
  sizes <- lengths(groups)
  log_nperm <- sum(lfactorial(sizes))
  max_log <- log(.Machine$double.xmax) - 1
  if(is.na(log_nperm) || log_nperm > max_log)
    nperm <- Inf
  else
    nperm <- round(exp(log_nperm))
  list(groups = groups,
       n_permutations = nperm,
       multi_units = sum(sizes > 1))
}


.permute_within_blocks <- function(groups){
  permuted <- integer(sum(lengths(groups)))
  for(block in groups){
    if(length(block) == 1){
      permuted[block] <- block
    } else {
      permuted[block] <- block[sample(length(block))]
    }
  }
  permuted
}


.all_permutations <- function(vec){
  n <- length(vec)
  if(n <= 1)
    return(matrix(vec, nrow = 1))
  total <- factorial(n)
  res <- matrix(0, nrow = total, ncol = n)
  idx <- 1L
  permute_rec <- function(prefix, remaining){
    if(length(remaining) == 0){
      res[idx, ] <<- prefix
      idx <<- idx + 1L
      return()
    }
    for(j in seq_along(remaining))
      permute_rec(c(prefix, remaining[j]), remaining[-j])
  }
  permute_rec(integer(0), vec)
  res
}


.enumerate_permutation_indices <- function(groups, group_perms, callback){
  order_template <- seq_len(sum(lengths(groups)))
  dims <- vapply(group_perms, nrow, integer(1))
  if(any(dims == 0))
    return()
  grid <- expand.grid(lapply(dims, seq_len), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  combos <- as.matrix(grid)
  for(row_idx in seq_len(nrow(combos))){
    order <- order_template
    for(group_idx in seq_along(groups)){
      ids <- groups[[group_idx]]
      order[ids] <- group_perms[[group_idx]][combos[row_idx, group_idx], ]
    }
    callback(order)
  }
}


.format_number <- function(x){
  if(!is.finite(x))
    return("a very large number")
  prettyNum(x, big.mark = ",", scientific = FALSE)
}


.run_cpp_permutation_with_progress <- function(DD, DR, n_perm, format){
  # Call C++ permutation kernel in chunks with progress bar
  vals <- numeric(n_perm)
  pb <- progress::progress_bar$new(total = n_perm, format = format)
  chunk_size <- min(250L, n_perm)
  from <- 1L
  while(from <= n_perm){
    n_chunk <- min(chunk_size, n_perm - from + 1L)
    to <- from + n_chunk - 1L
    # For approximate permutation, use a single global block
    groups_global <- list(seq_len(nrow(DR)))
    vals[from:to] <- .permutation_ssq_kernel_cpp(DD, DR, groups_global, n_chunk)
    pb$tick(n_chunk)
    from <- to + 1L
  }
  vals
}