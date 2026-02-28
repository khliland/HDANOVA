#' Rotation test for HDANOVA
#'
#' @description Rotation testing for HDANOVA. This function performs
#' random orthogonal rotations of exchangeable residual units for each approved
#' effect and adds the resulting null distributions to the \code{hdanova} object.
#'
#' @param object A \code{hdanova} object.
#' @param rotate Number of random rotations to perform (default = 1000).
#' @param unique.digits Number of digits used when rounding rotation SSQ values before
#' checking uniqueness (default = 12). Set to \code{NULL} to disable this warning.
#' @param unique.frac Minimum fraction of unique rounded SSQ values required to avoid
#' warning (default = 0.95). Set to \code{NULL} to disable this warning.
#' @param block.type Rotation blocking strategy. \code{"denominator"} (default)
#' rotates within denominator-compatible exchangeable blocks. \code{"global"}
#' rotates across all observations.
#'
#' @returns An updated \code{hdanova} object with rotation-test results stored in
#' \code{object$permute} for compatibility with existing summary and plotting tools.
#'
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic HDANOVA model with two factors
#' mod <- hdanova(assessment ~ candy + assessor, data=candies)
#'
#' # Rotation test
#' modRot <- rotation(mod)
#' summary(modRot)
#'
#' # Plot null distribution for "candy" effect
#' rotationplot(modRot, factor="candy")
#'
#' @export
rotation <- function(object,
                     rotate=1000,
                     unique.digits=12,
                     unique.frac=0.95,
                     block.type=c("denominator", "global")){
  block.type <- match.arg(block.type)
  check_unique_warning <- !is.null(unique.digits) && !is.null(unique.frac)
  if(check_unique_warning){
    if(!is.numeric(unique.digits) || length(unique.digits) != 1 || is.na(unique.digits))
      stop("'unique.digits' must be a single numeric value or NULL.")
    if(!is.numeric(unique.frac) || length(unique.frac) != 1 || is.na(unique.frac) ||
       unique.frac < 0 || unique.frac > 1)
      stop("'unique.frac' must be a single numeric value in [0, 1] or NULL.")
    unique.digits <- as.integer(unique.digits)
  }

  ssqa <- pvals <- numeric(0)
  ssqaperm <- list()
  rot_warnings <- character(0)
  exchangeable_units <- list()
  rotation_orders <- list()

  if(!(is.logical(rotate) && !rotate)){
    if(is.logical(rotate)){
      rotate <- 1000
      if(interactive())
        message("Defaulting to 1000 rotations\n")
    }
    if(!is.numeric(rotate) || length(rotate) != 1 || is.na(rotate) || rotate < 1)
      stop("'rotate' must be a positive integer.")
    rotate <- as.integer(rotate)

    approved <- object$more$approved
    eff_names <- object$more$effs[approved]
    ssqa <- pvals <- numeric(length(approved))
    ssqaperm <- vector("list", length(approved))
    exchangeable_units <- vector("list", length(approved))
    rotation_orders <- vector("list", length(approved))
    names(ssqa) <- names(pvals) <- names(ssqaperm) <-
      names(exchangeable_units) <- names(rotation_orders) <- eff_names
    cpp_available <- exists(".rotation_ssq_kernel_cpp", mode = "function")
    cpp_failed <- FALSE

    for(i in seq_along(approved)){
      a <- object$more$approvedAB[i]
      eff_name <- object$more$effs[a]
      D <- object$X[, object$more$assign %in% object$more$approvedComb[[names(a)]], drop=FALSE]
      DD <- D %*% pracma::pinv(D)
      DR <- object$more$LS_aug[[eff_name]]
      ssqa[eff_name] <- norm(DD %*% DR, "F")^2

      if(block.type == "global"){
        groups <- list(seq_len(object$more$N))
      } else {
        effect_vars <- .rotation_effect_variable_names(object, eff_name)
        block_info <- .rotation_build_blocks(object, eff_name, effect_vars)
        groups <- block_info$groups
      }

      sizes <- lengths(groups)
      active_groups <- groups[sizes > 1]
      exchangeable_units[[eff_name]] <- list(
        total_blocks = length(groups),
        exchangeable_blocks = sum(sizes > 1),
        block_sizes = sizes
      )

      if(all(sizes <= 1)){
        rot_warnings <- c(rot_warnings,
                          paste0("Effect '", eff_name, "' has no exchangeable units; rotation test skipped."))
        ssqaperm[[eff_name]] <- numeric(0)
        rotation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
        pvals[eff_name] <- NA
        next
      }

      vals <- NULL
      if(cpp_available && !cpp_failed){
        cpp_res <- tryCatch(
          list(values = .run_cpp_rotation_with_progress(
            DD = DD,
            DR = DR,
            groups = active_groups,
            n_rot = rotate,
            format = paste0("  Rotating ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)")
          ),
               error = NULL),
          error = function(e)
            list(values = NULL, error = e)
        )
        vals <- cpp_res$values
        if(!is.null(cpp_res$error)){
          cpp_failed <- TRUE
          rot_warnings <- c(rot_warnings,
                            paste0("C++ rotation kernel failed for effect '", eff_name,
                                   "'; falling back to R implementation. ", cpp_res$error$message))
        }
      }
      if(is.null(vals)){
        vals <- numeric(rotate)
        pb <- progress_bar$new(total = rotate,
                               format = paste0("  Rotating ", eff_name, " (", i, "/", length(approved), ") [:bar] :percent (:eta)"))
        for(rot in seq_len(rotate)){
          DR_rot <- .rotate_within_blocks(DR, active_groups)
          vals[rot] <- norm(DD %*% DR_rot, "F")^2
          pb$tick()
        }
      }

      ssqaperm[[eff_name]] <- vals
      rotation_orders[[eff_name]] <- matrix(integer(0), nrow = 0, ncol = object$more$N)
      if(check_unique_warning && length(vals) > 0){
        unique_ssq <- length(unique(round(vals, unique.digits)))
        if(unique_ssq < length(vals) * unique.frac){
          rot_warnings <- c(rot_warnings,
                            paste0("Effect '", eff_name, "' produced only ",
                                   .rotation_format_number(unique_ssq),
                                   " unique SSQ values (rounded to ", unique.digits,
                                   " decimals) out of ",
                                   .rotation_format_number(length(vals)),
                                   " evaluated rotations."))
        }
      }
      pvals[eff_name] <- sum(vals > ssqa[eff_name])/(length(vals))
    }
  }

  if(length(rot_warnings) > 0)
    warning(paste(unique(rot_warnings), collapse = "\n"))

  object$permute <- list(ssqa=ssqa,
                         ssqaperm=ssqaperm,
                         pvalues=pvals,
                         permutations=rotate,
                         exchangeable=exchangeable_units,
                         orders=rotation_orders,
                         method="rotation",
                         block.type=block.type)
  object
}

.rotate_within_blocks <- function(mat, groups){
  for(block in groups){
    m <- length(block)
    if(m <= 1)
      next
    Q <- .random_orthogonal_matrix(m)
    mat[block, ] <- Q %*% mat[block, , drop=FALSE]
  }
  mat
}

.random_orthogonal_matrix <- function(n){
  z <- matrix(stats::rnorm(n * n), nrow = n)
  qr_z <- qr(z)
  q <- qr.Q(qr_z)
  d <- sign(diag(qr_z$qr))
  d[d == 0] <- 1
  sweep(q, 2, d, `*`)
}

.rotation_effect_variable_names <- function(object, eff_name){
  combo <- object$more$approvedComb[[eff_name]]
  if(is.null(combo))
    combo <- which(object$more$effs == eff_name)
  if(length(combo) == 0)
    combo <- which(object$more$effs %in% eff_name)
  names <- object$more$effs[combo]
  unique(unlist(lapply(names, .rotation_parse_effect_terms)))
}

.rotation_parse_effect_terms <- function(name){
  if(is.null(name))
    return(character(0))
  clean <- gsub("\\s+", "", name)
  if(clean == "" || clean == "(Intercept)")
    return(character(0))
  terms <- strsplit(clean, "[:+\\|]", perl = TRUE)[[1]]
  terms[terms != ""]
}

.rotation_build_blocks <- function(object, eff_name, effect_vars){
  mf <- object$model.frame
  if(ncol(mf) > 0)
    mf <- mf[setdiff(names(mf), names(mf)[1])]
  candidate <- intersect(names(mf), object$more$effs)
  candidate <- setdiff(candidate, eff_name)
  candidate <- candidate[!candidate %in% "(Intercept)"]
  eligible <- vapply(candidate, function(name){
    vars <- .rotation_parse_effect_terms(name)
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
  list(groups = groups)
}

.rotation_format_number <- function(x){
  if(!is.finite(x))
    return("a very large number")
  prettyNum(x, big.mark = ",", scientific = FALSE)
}

.run_cpp_rotation_with_progress <- function(DD, DR, groups, n_rot, format){
  vals <- numeric(n_rot)
  pb <- progress_bar$new(total = n_rot, format = format)
  chunk_size <- min(250L, n_rot)
  from <- 1L
  while(from <= n_rot){
    n_chunk <- min(chunk_size, n_rot - from + 1L)
    to <- from + n_chunk - 1L
    vals[from:to] <- .rotation_ssq_kernel_cpp(DD, DR, groups, n_chunk)
    pb$tick(n_chunk)
    from <- to + 1L
  }
  vals
}
