#' ANOVA Principal Component Analysis - APCA
#'
#' @description APCA function for fitting ANOVA Principal Component Analysis models.
#'
#' @param formula Model formula accepting a single response (block) and predictors.
#' @param data The data set to analyse.
#' @param ... Additional parameters for the asca_fit function.
#'
#' @return An object of class \code{apca}, inheriting from the general \code{asca} class.
#' TODO: Link to ASCA documentation, add reference
#' @export
#'
#' @examples
#' data(candies)
#' ap <- apca(assessment ~ candy, data=candies)
#' scoreplot(ap)
apca <- function(formula, data, ...){
  object <- asca_fit(formula, data, add_error = TRUE, ...)
  object$call <- match.call()
  class(object) <- c("apca", class(object))
  object
}
