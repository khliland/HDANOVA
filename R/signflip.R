#' @title Flip signs of a component/factor combination in a SCA/PCA object
#'
#' @description This function flips the sign of a selected component in a selected factor of an \code{asca} object.
#' This affects both scores, loadings and projected data.
#'
#' @param object \code{asca} object.
#' @param factor \code{integer/character} for selecting a model factor.
#' @param comp \code{integer} for selected component.
#'
#' @returns
#' An \code{asca} object with the sign of the selected component flipped.
#'
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic HDANOVA model with two factors
#' mod <- hdanova(assessment ~ candy + assessor, data=candies)
#' mod <- sca(mod)
#' old.par <- par(mfrow=c(1,2), mar=c(4,4,1,1))
#' scoreplot(mod, factor="candy")
#' loadingplot(mod, factor="candy")
#' par(old.par)
#'
#' # Flip the sign of the first component of the candy factor
#' mod <- signflip(mod, factor="candy", comp=1)
#' old.par <- par(mfrow=c(1,2), mar=c(4,4,1,1))
#' scoreplot(mod, factor="candy")
#' loadingplot(mod, factor="candy")
#' par(old.par)
#' @export
signflip <- function(object, factor, comp){
  if(length(factor) != 1 || length(comp) != 1)
    stop("Use this function for one factor and component at a time")
  if(is.null(object$loadings) || is.null(object$scores) || is.null(object$projected)){
    stop("Object contains no SCA/PCA results.")
  }
  object$loadings[[factor]][,comp] <- -object$loadings[[factor]][,comp]
  object$scores[[factor]][,comp] <- -object$scores[[factor]][,comp]
  object$projected[[factor]][,comp] <- -object$projected[[factor]][,comp]
  object
}
