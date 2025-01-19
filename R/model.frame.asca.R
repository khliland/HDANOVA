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
#' @returns
#' @export
#'
#' @examples
model.frame.asca <- function(formula, ...){
  return(formula$model)
}

#' @rname model.frame.asca
#' @export
model.matrix.asca <- function(object, ...){
  return(object$X)
}
