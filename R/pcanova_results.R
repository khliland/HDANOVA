#' @title PC-ANOVA Result Methods
#' @name pcanova_results
#' @aliases pcanova_results print.pcanova summary.pcanova projections projections.pcanova print.summary.pcanova loadings.pcanova scores.pcanova
#'
#' @description Standard result computation and extraction functions for ASCA (\code{\link{pcanova}}).
#'
#' @details Usage of the functions are shown using generics in the examples in \code{\link{pcanova}}.
#' Explained variances are available (block-wise and global) through \code{blockexpl} and \code{print.rosaexpl}.
#' Object printing and summary are available through:
#' \code{print.pcanova} and \code{summary.pcanova}.
#' Scores and loadings have their own extensions of \code{scores()} and \code{loadings()} through
#' \code{scores.pcanova} and \code{loadings.pcanova}. Special to ASCA is that scores are on a
#' factor level basis, while back-projected samples have their own function in \code{projections.pcanova}.
#'
#' @param object \code{pcanova} object.
#' @param x \code{pcanova} object.
#' @param factor \code{integer/character} for selecting a model factor.
#' @param digits \code{integer} number of digits for printing.
#' @param ... additional arguments to underlying methods.
#'
#' @return Returns depend on method used, e.g. \code{projections.pcanova} returns projected samples,
#' \code{scores.pcanova} return scores, while print and summary methods return the object invisibly.
#'
#' @references Luciano G, Næs T. Interpreting sensory data by combining principal
#' component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167-175.
#'
#' @seealso TODO Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{pcanova}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for plotting are found in \code{\link{pcanova_plots}}.
#'
#' @export
summary.pcanova <- function(object, ...){
  anos <- object$anovas
  x <- list(dat=anos, fit.type=object$fit.type, explvar=object$explvar)
  class(x) <- c('summary.pcanova')
  x
}

#' @rdname pcanova_results
#' @export
print.summary.pcanova <- function(x, digits=2, ...){
  anos <- x$dat
  names(anos) <- paste0(names(anos), ' (', round(x$explvar, digits), '%)')
  cat("PC-ANOVA fitted using", x$fit.type, "\n")
  print(anos)
  invisible(anos)
}
