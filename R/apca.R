#' ANOVA Principal Component Analysis - APCA
#'
#' @description APCA function for fitting ANOVA Principal Component Analysis models.
#'
#' @param formula Model formula accepting a single response (block) and predictors.
#' @param data The data set to analyse.
#' @param add_error Add error to LS means (default = TRUE).
#' @param ... Additional parameters for the asca_fit function.
#'
#' @return An object of class \code{apca}, inheriting from the general \code{asca} class.
#' Further arguments and plots can be found in the \code{\link{asca}} documentation.
#' @references
#' Harrington, P.d.B., Vieira, N.E., Espinoza, J., Nien, J.K., Romero, R., and Yergey, A.L. (2005)
#' Analysis of variance–principal component analysis: A soft tool for proteomic discovery.
#' Analytica chimica acta, 544 (1-2), 118–127.
#' @export
#'
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{asca_fit}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#' @examples
#' data(candies)
#' ap <- apca(assessment ~ candy, data=candies)
#' scoreplot(ap)
#'
#' # Numeric effects
#' candies$num <- eff <- 1:165
#' mod <- apca(assessment ~ candy + assessor + num, data=candies)
#' summary(mod)
#' scoreplot(mod, factor=3, gr.col=rgb(eff/max(eff), 1-eff/max(eff),0), pch.scores="x")
apca <- function(formula, data, add_error = TRUE, ...){
  object <- asca_fit(formula, data, add_error = add_error, ...)
  object$call <- match.call()
  class(object) <- c("apca", class(object))
  object
}
