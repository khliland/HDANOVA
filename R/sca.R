#' Simultaneous Component Analysis
#'
#' @description This function performs Simultaneous Component Analysis (SCA) on a \code{hdanova} object.
#'
#' @param object A \code{hdanova} object.
#'
#' @returns An updated \code{hdanova} object with SCA results.
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic HDANOVA model with two factors
#' mod <- hdanova(assessment ~ candy + assessor, data=candies)
#' mod <- sca(mod)
#' scoreplot(mod)
#'
#' @export
sca <- function(object){
  ########################## SCA ##########################
  # SCAs
  scores <- loadings <- projected <- singulars <- list()
  for(i in object$more$approved){
    maxDiri <- min(Rank(object$LS[[object$more$effs[i]]]),object$more$maxDir[i])
    if(object$more$pca.in != 0)
      maxDiri <- min(maxDiri, object$more$pca.in)
    if(object$add_error)
      maxDiri <- min(object$more$N-1, object$more$p)
    if(maxDiri == 0)
      stop(paste0("Effect '", object$more$effs[i], "' has no estimable levels"))
    pcai <- .pca(object$LS[[object$more$effs[i]]], ncomp=maxDiri, proj=object$error[[object$more$effs[i]]])
    scores[[object$more$effs[i]]] <- pcai$scores
    loadings[[object$more$effs[i]]] <- pcai$loadings
    projected[[object$more$effs[i]]] <- pcai$projected
    singulars[[object$more$effs[i]]] <- pcai$singulars

    if(object$more$pca.in!=0){ # Transform back if PCA on Y has been performed
      loadings[[object$more$effs[i]]] <- object$Ypca$pca$loadings[,1:object$more$pca.in,drop=FALSE] %*% loadings[[object$more$effs[i]]]
      dimnames(loadings[[object$more$effs[i]]]) <- list(colnames(object$Y), paste("Comp", 1:maxDiri, sep=" "))
    }
  }
  # SCA of residuals
  maxDirRes <- min(object$more$N-1,object$more$p)
  if(object$more$pca.in != 0)
    maxDirRes <- min(maxDirRes, object$more$pca.in)
  pcaRes <- .pca(object$residuals, ncomp=maxDirRes)
  scores[["Residuals"]] <- pcaRes$scores
  loadings[["Residuals"]] <- pcaRes$loadings
  projected[["Residuals"]] <- pcaRes$projected
  singulars[["Residuals"]] <- pcaRes$singulars

  ########################## Return ##########################
  object$scores <- scores
  object$loadings <- loadings
  object$projected <- projected
  object$singulars <- singulars
  class(object) <- c("asca", class(object))
  return(object)
}
