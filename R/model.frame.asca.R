#' @name model.frame.asca
#' @aliases model.frame.asca model.matrix.asca
#' @title Model Frame and Model Matrix for ASCA-like Models
#'
#' @description
#' Extraction functions to retrieve the \code{model.frame} and \code{model.matrix} of an \code{asca} object.
#'
#' @param formula The \code{asca} object.
#' @param object The \code{asca} object.
#' @param ... Not implemented
#'
#' @returns A \code{data.frame} or \code{matrix} object.
#' @export
#'
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic ASCA model with two factors
#' mod <- asca(assessment ~ candy + assessor, data=candies)
#'
#' # Extract model frame and model matrix
#' mf <- model.frame(mod)
#' head(mf)
#' mm <- model.matrix(mod)
#' par.old <- par(mar=c(3,3,3,1), mgp=c(1,0.7,0))
#' image(t(mm[seq(165,1,-1),]), main="Model Matrix", xlab="dummy values", ylab="samples",
#'      axes=FALSE)
#' par(par.old)
model.frame.asca <- function(formula, ...){
  return(formula$model)
}

#' @rdname model.frame.asca
#' @export
model.matrix.asca <- function(object, ...){
  return(object$X)
}
