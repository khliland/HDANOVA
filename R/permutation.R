#' Permutation for HDANOVA
#'
#' @description Permutation testing for HDANOVA. This function performes
#' permutation testing for the effects in the HDANOVA model and adds them to the
#' \code{hdanova} object.
#'
#' @param object A \code{hdanova} object.
#' @param permute Number of permutations to perform (default = 1000).
#' @param perm.type Type of permutation to perform, either "approximate" or "exact" (default = "approximate").
#'
#' @returns An updated \code{hdanova} object with permutation results.
#'
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic HDANOVA model with two factors
#' mod <- hdanova(assessment ~ candy + assessor, data=candies)
#' mod <- permutation(mod)
#' summary(mod)
#'
#' @export
permutation <- function(object,
                        permute=1000,
                        perm.type=c("approximate","exact")){
  ########################## Permutation ##########################
  # Permutation testing
  if(!(is.logical(permute) && !permute)){
    # Default to 1000 permutations             ------------- Flytt testing til ASCA og venner
    if(is.logical(permute)){
      permute <- 1000
      if(interactive())
        message("Defaulting to 1000 permutations\n")
    }
    if(perm.type[1] == "approximate"){
      ssqa <- pvals <- numeric(length(object$more$approved))
      ssqaperm <- lapply(1:length(object$more$approved),function(i)"NA")
      names(ssqa) <- names(ssqaperm) <- names(pvals) <- object$more$effs[object$more$approved]
      for(i in 1:length(object$more$approved)){
        a <- object$more$approvedAB[i]
        perms <- numeric(permute)
        # Subset of design matrix for effect a
        D <- object$X[, object$more$assign%in%object$more$approvedComb[[names(a)]], drop=FALSE]
        DD <- D %*% pracma::pinv(D)
        # Base ssq
        ssqa[object$more$effs[a]] <- norm(DD %*% object$more$LS_aug[[object$more$effs[a]]], "F")^2
        pb <- progress_bar$new(total = permute, format = paste0("  Permuting ", object$more$effs[a], " (", i,"/",length(object$more$approved),") [:bar] :percent (:eta)"))
        # Permuted ssqs
        for(perm in 1:permute){
          perms[perm] <- norm(DD %*% object$more$LS_aug[[object$more$effs[a]]][sample(object$more$N),], "F")^2
          pb$tick()
        }
        ssqaperm[[object$more$effs[a]]] <- perms
        pvals[object$more$effs[a]] <- sum(perms > ssqa[object$more$effs[a]])/(permute)
      }
    }
    if(perm.type[1] == "exact"){
      # lme4 vs LS
      if(!object$more$lme4 && !is.logical(object$more$REML)){
        stop("Exact permutations in lme4 models not implemented yet!")
      }
      # Loop over effects
      for(i in 1:length(object$more$approved)){
        stop("Exact permutations not implemented yet!")
        a <- object$more$approvedAB[i]
        perms <- numeric(permute)
        # Subset of design matrix for effect a
        #D <- M[, assign%in%approvedComb[[names(a)]], drop=FALSE]
        #DD <- D %*% pracma::pinv(D)
        # Find permissible units
        # Check balance, warning
        # Permute
      }
    }
  }
  ########################## Return ##########################
  object$permute <- list(ssqa=ssqa, ssqaperm=ssqaperm, pvalues=pvals, permutations=permute)
  return(object)
}
