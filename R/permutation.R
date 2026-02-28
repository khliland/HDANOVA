#' Permutation for HDANOVA
#'
#' @description Permutation testing for HDANOVA. This function performes
#' permutation testing for the effects in the HDANOVA model and adds them to the
#' \code{hdanova} object.
#'
#' @param object A \code{hdanova} object.
#' @param permute Number of permutations to perform (default = 10000).
#' @param perm.type Type of permutation to perform, either "approximate" or "exact" (default = "approximate").
#' @param unique.digits Number of digits used when rounding permutation SSQ values before
#' checking uniqueness (default = 12). Set to \code{NULL} to disable this warning.
#' @param unique.frac Minimum fraction of unique rounded SSQ values required to avoid
#' warning (default = 0.95). Set to \code{NULL} to disable this warning.
#' @param exhaustive.warn Logical; if \code{TRUE} (default), warn when exact permutation
#' uses exhaustive enumeration with fewer permutations than requested.
#'
#' @returns An updated \code{hdanova} object with permutation results.
#'
#' @examples
#' # Load candies data
#' data(candies)
#'
#' ## Basic HDANOVA model with two factors
#' mod <- hdanova(assessment ~ candy + assessor, data=candies)
#'
#' ## Approximate permutation
#' modApprox <- permutation(mod)
#' summary(modApprox)
#'
#' # Plot permutation distribution for "candy" effect
#' permutationplot(modApprox, factor="candy")
#'
#' ## Exact permutation (warning if too few exchangeable units)
#' modExact  <- permutation(mod, perm.type="exact")
#' summary(modExact)
#'
#' # Reduced candy data (first two levels of each effect, two replicates)
#' reduced_candies <- candies[candies$candy %in% levels(candies$candy)[1:2] &
#'      candies$assessor %in% levels(candies$assessor)[1:2], ][-seq(1,12,by=3),]
#' mod_reduced <- hdanova(assessment ~ candy + assessor, data=reduced_candies)
#' mod_reducedApprox <- permutation(mod_reduced)
#' mod_reducedExact  <- permutation(mod_reduced, perm.type="exact")
#' par.old <- par(mfrow=c(2,1))
#' permutationplot(mod_reducedApprox, factor="assessor", main="Approximate permutation")
#' permutationplot(mod_reducedExact, factor="assessor", main="Exact permutation")
#' par(par.old)
#'
#' # Check how many exchangeable units were available (minimum per effect)
#' lapply(modExact$permute$exchangeable, function(x)min(x$block_sizes))
#'
#' # Caldana data (combined effects and exact permutation)
#' data(caldana)
#' mod.comb <- asca(compounds ~ time + comb(light + time:light), data=caldana)
#' mod.comb <- permutation(mod.comb, perm.type="exact")
#' summary(mod.comb)
#'
#' @export
permutation <- function(object,
                        permute=1000,
                        perm.type=c("approximate","exact"),
                        unique.digits=12,
                        unique.frac=0.95,
                        exhaustive.warn=TRUE){
  ########################## Permutation ##########################
    perm.type <- match.arg(perm.type)
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
    cpp_available <- exists(".permutation_ssq_kernel_cpp", mode = "function")
    cpp_failed <- FALSE

    if(perm.type == "approximate"){
      for(i in seq_along(approved)){
        a <- object$more$approvedAB[i]
        eff_name <- object$more$effs[a]
        perms <- NULL
        D <- object$X[, object$more$assign %in% object$more$approvedComb[[names(a)]], drop=FALSE]
        DD <- D %*% pracma::pinv(D)
        DR <- object$more$LS_aug[[eff_name]]
        ssqa[eff_name] <- norm(DD %*% DR, "F")^2
        if(cpp_available && !cpp_failed){
          cpp_res <- tryCatch(
            list(values = .run_cpp_permutation_with_progress(
              DD = DD,
              DR = DR,
              groups = list(seq_len(object$more$N)),
              n_perm = permute,
              format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)")
            ),
                 error = NULL),
            error = function(e)
              list(values = NULL, error = e)
          )
          if(!is.null(cpp_res$error)){
            cpp_failed <- TRUE
            perm_warnings <- c(perm_warnings,
                               paste0("C++ permutation kernel failed for effect '", eff_name,
                                      "'; falling back to R implementation. ", cpp_res$error$message))
          } else {
            perms <- cpp_res$values
          }
        }
        if(is.null(perms) || length(perms) != permute){
          perms <- numeric(permute)
          pb <- progress_bar$new(total = permute,
                                 format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)"))
          for(perm in seq_len(permute)){
            perms[perm] <- norm(DD %*% DR[sample(object$more$N),], "F")^2
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
        pvals[eff_name] <- sum(perms > ssqa[eff_name])/(permute)
      }
    } else {
      for(i in seq_along(approved)){
        a <- object$more$approvedAB[i]
        eff_name <- object$more$effs[a]
        D <- object$X[, object$more$assign %in% object$more$approvedComb[[names(a)]], drop=FALSE]
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
          pvals[eff_name] <- NA
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
            perms[counter] <<- norm(DD %*% DR[order,], "F")^2
          })
          permutation_orders[[eff_name]] <- order_matrix
        } else {
          perms <- NULL
          if(cpp_available && !cpp_failed){
            cpp_res <- tryCatch(
              list(values = .run_cpp_permutation_with_progress(
                DD = DD,
                DR = DR,
                groups = block_info$groups,
                n_perm = permute,
                format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)")
              ),
                   error = NULL),
              error = function(e)
                list(values = NULL, error = e)
            )
            if(!is.null(cpp_res$error)){
              cpp_failed <- TRUE
              perm_warnings <- c(perm_warnings,
                                 paste0("C++ permutation kernel failed for effect '", eff_name,
                                        "'; falling back to R implementation. ", cpp_res$error$message))
            } else {
              perms <- cpp_res$values
            }
          }
          if(is.null(perms) || length(perms) != permute){
            perms <- numeric(permute)
            pb <- progress_bar$new(total = permute,
                                   format = paste0("  Permuting ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)"))
            for(perm in seq_len(permute)){
              idx <- .permute_within_blocks(block_info$groups)
              perms[perm] <- norm(DD %*% DR[idx,], "F")^2
              pb$tick()
            }
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
        pvals[eff_name] <- sum(perms > ssqa[eff_name])/(length(perms))
      }
    }
  }

  if(length(perm_warnings) > 0)
    warning(paste(unique(perm_warnings), collapse = "\n"))

  ########################## Return ##########################
  object$permute <- list(ssqa=ssqa,
                         ssqaperm=ssqaperm,
                         pvalues=pvals,
                         permutations=permute,
                         exchangeable=exchangeable_units,
                         orders=permutation_orders,
                         method="permutation",
                         perm.type=perm.type)
  return(object)
}

#################### Helper utilities ####################
.effect_variable_names <- function(object, eff_name){
  # Extract the categorical/interaction names that compose the effect;
  # denominator MS for this effect is determined by everything not orthogonal
  # to these terms, so we want to avoid permuting within any variable that
  # contributes variability to the denominator (the exchangeable units).
  combo <- object$more$approvedComb[[eff_name]]
  if(is.null(combo))
    combo <- which(object$more$effs == eff_name)
  if(length(combo) == 0)
    combo <- which(object$more$effs %in% eff_name)
  names <- object$more$effs[combo]
  unique(unlist(lapply(names, .parse_effect_terms)))
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
  # Determine exchangeable units by splitting observations along variables that
  # do not contribute to the denominator MS of the tested effect. The remaining
  # variables define blocks (strata) where residual variation comes solely from
  # the error term, matching the notion that the denominator MS identifies the
  # components that can be permuted.
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
    keys <- do.call(paste, c(data, sep="\r"))
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
    for(j in seq_along(remaining)){
      permute_rec(c(prefix, remaining[j]), remaining[-j])
    }
  }
  permute_rec(integer(0), vec)
  res
}

.enumerate_permutation_indices <- function(groups, group_perms, callback){
  order_template <- seq_len(sum(lengths(groups)))
  dims <- vapply(group_perms, nrow, integer(1))
  if(any(dims == 0))
    return()
  # Cartesian product of permutation indices across blocks replaces recursion.
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

.run_cpp_permutation_with_progress <- function(DD, DR, groups, n_perm, format){
  vals <- numeric(n_perm)
  pb <- progress_bar$new(total = n_perm, format = format)
  chunk_size <- min(250L, n_perm)
  from <- 1L
  while(from <= n_perm){
    n_chunk <- min(chunk_size, n_perm - from + 1L)
    to <- from + n_chunk - 1L
    vals[from:to] <- .permutation_ssq_kernel_cpp(DD, DR, groups, n_chunk)
    pb$tick(n_chunk)
    from <- to + 1L
  }
  vals
}
