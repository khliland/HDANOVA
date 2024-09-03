#' @title PC-ANOVA Result Methods
#' @name pcanova_plots
#' @aliases pcanova_plots scoreplot.pcanova loadingplot.pcanova
#'
#' @description Various plotting procedures for \code{\link{pcanova}} objects.
#'
#' @details Usage of the functions are shown using generics in the examples in \code{\link{pcanova}}.
#' Plot routines are available as
#' \code{scoreplot.pcanova} and \code{loadingplot.pcanova}.
#'
#' @param object \code{pcanova} object.
#' @param factor \code{integer/character} for selecting a model factor.
#' @param comps \code{integer} vector of selected components.
#' @param ... additional arguments to underlying methods.
#'
#' @return The plotting routines have no return.
#'
#' @references Luciano G, Næs T. Interpreting sensory data by combining principal
#' component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167-175.
#'
#' @seealso TODO Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{pcanova}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results are found in \code{\link{pcanova_results}}.
#'
#' @export
scoreplot.pcanova <- function(object, factor = 1, comps = 1:2, col = "factor", ...){
  # Remove too high component numbers
  comps <- comps[comps <= length(object$anovas)]
  if(col == "factor"){
    levels <- model.frame(object$models[[1]])[,-1,drop=FALSE]
    col <- levels[[factor]]
  }
  scoreplot(object$scores, comps=comps, col=col, ...)
}
