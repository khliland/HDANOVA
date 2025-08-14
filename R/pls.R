#' Partial Least Squares (PLS) for HDANOVA
#'
#' @description This function performs Partial Least Squares (PLS) on a \code{hdanova}.
#'
#' @param object A \code{hdanova} object.
#'
#' @returns An updated \code{hdanova} object with PLS results.
#' @details For residuals, PCA is performed instead of PLS as there is no natural response.
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic HDANOVA model with two factors
#' mod <- hdanova(assessment ~ candy + assessor, data=candies)
#' mod <- pls(mod)
#' scoreplot(mod)
#'
#' @export
pls <- function(object){
  scores <- loadings <- projected <- list()
  for(i in object$more$approved){
    maxDiri <- min(Rank(object$LS[[object$more$effs[i]]]),object$more$maxDir[i])
    if(object$more$pca.in != 0)
      maxDiri <- min(maxDiri, object$more$pca.in)
    if(object$add_error)
      maxDiri <- min(object$more$N-1, object$more$p)
    if(maxDiri == 0)
      stop(paste0("Effect '", object$more$effs[i], "' has no estimable levels"))
    # Check for combined effect
    if(object$eff_combined[names(which(object$more$approved==i))]){
      combs <- object$more$approvedComb[[names(which(object$more$approved==i))]]
      if(inherits(object$model.frame[[object$more$effs[i]]], "numeric") ||
         inherits(object$model.frame[[object$more$effs[i]]], "integer"))
        Y <- object$model.frame[[combs[1]]]
      else
        Y <- dummycode(object$model.frame[[combs[1]]], contrast="contr.treatment", drop = FALSE)
      for(j in 2:length(combs)){
        if(inherits(object$model.frame[[object$more$effs[i]]], "numeric") ||
           inherits(object$model.frame[[object$more$effs[i]]], "integer"))
          Y <- cbind(Y, object$model.frame[[combs[j]]])
        else
          Y <- cbind(Y, dummycode(object$model.frame[[combs[j]]], contrast="contr.treatment", drop = FALSE))
      }
    } else {
      if(inherits(object$model.frame[[object$more$effs[i]]], "numeric") ||
         inherits(object$model.frame[[object$more$effs[i]]], "integer"))
        Y <- object$model.frame[[object$more$effs[i]]]
      else
        Y <- dummycode(object$model.frame[[object$more$effs[[i]]]], contrast="contr.treatment", drop = FALSE)
    }
    dat <- data.frame(X = I(object$LS[[object$more$effs[i]]]),
                      Y = I(Y))
    plsi <- plsr(Y~X, data=dat, ncomp=maxDiri)
#    plsi <- .pca(object$LS[[object$more$effs[i]]], ncomp=maxDiri, proj=object$error[[object$more$effs[i]]])
    scores[[object$more$effs[i]]] <- plsi$scores
    attr(scores[[object$more$effs[i]]], "explvar") <- pls::explvar(plsi)
    loadings[[object$more$effs[i]]] <- plsi$loadings
    projected[[object$more$effs[i]]] <- object$residuals%*%plsi$projection + plsi$scores

    if(object$more$pca.in!=0){ # Transform back if PCA on Y has been performed
      loadings[[object$more$effs[i]]] <- object$Ypca$pca$loadings[,1:object$more$pca.in,drop=FALSE] %*% loadings[[object$more$effs[i]]]
      dimnames(loadings[[object$more$effs[i]]]) <- list(colnames(object$Y), paste("Comp", 1:maxDiri, sep=" "))
    }
  }
  # PCA of residuals
  maxDirRes <- min(object$more$N-1,object$more$p)
  if(object$more$pca.in != 0)
    maxDirRes <- min(maxDirRes, object$more$pca.in)
  pcaRes <- .pca(object$residuals, ncomp=maxDirRes)
  scores[["Residuals"]] <- pcaRes$scores
  loadings[["Residuals"]] <- pcaRes$loadings
  projected[["Residuals"]] <- pcaRes$projected
#  singulars[["Residuals"]] <- pcaRes$singulars

  ########################## Return ##########################
  object$scores <- scores
  object$loadings <- loadings
  object$projected <- projected
#  object$singulars <- singulars
  if(inherits(object,"asca"))
    class(object) <- c("apls", class(object))
  else
    class(object) <- c("apls","asca", class(object))
  return(object)
}
