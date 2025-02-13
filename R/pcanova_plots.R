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
#' @param col \code{character} for selecting a factor to use for colouring
#' (default = first factor) or ordinary colour specifications.
#' @param ... additional arguments to underlying methods.
#'
#' @return The plotting routines have no return.
#'
#' @references Luciano G, Næs T. Interpreting sensory data by combining principal
#' component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167-175.
#'
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{hdanova}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
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
