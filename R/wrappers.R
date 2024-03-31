#' Principal Response Curves
#' @description
#' Wrapper for the \code{\link[vegan]{prc}} function to allow for formula input.
#'
#' @param formula Model formula accepting a single response (block) and predictors.
#' If no predictor is called 'time', time is assumed to be the second predictor.
#' @param data The data set to analyse.
#' @param ... Additional arguments to \code{\link[vegan]{prc}}.
#'
#' @return An object of class \code{prc}.
#' @export
#'
#' @examples
#' data(caldana)
#' (pr <- prc(compounds ~ light * time, caldana))
#' summary(pr)
prc <- function(formula, data, ...) {
  # Extract response and predictors from formula
  mf <- model.frame(formula,data)
  response <- model.response(mf)
  if(any(colnames(mf)=="time")){
    time      <- mf$time
    treatment <- mf[[setdiff(colnames(mf), "time")[2]]]
  } else {
    treatment <- mf[[2]]
    time      <- mf[[3]]
  }
  requireNamespace("vegan", quietly = TRUE)
  vegan::prc(response, treatment, time, ...)
}


#' Permutation Based MANOVA - PERMANOVA
#' @description
#' Wrapper for the \code{\link[vegan]{adonis2}} function to allow ordinary formula input.
#'
#' @param formula Model formula accepting a single response matrix and predictors.
#' See details in \code{\link[vegan]{adonis2}}.
#' @param data The data set to analyse.
#' @param ... Additional arguments to \code{\link[vegan]{adonis2}}.
#'
#' @return An ANOVA table with permutation-based p-values.
#' @export
#'
#' @examples
#' data(caldana)
#' (pr <- permanova(compounds ~ light * time, caldana))
permanova <- function(formula, data, ...){
  # Copy data to current environment
  ca <- match.call()
  assign(as.character(ca$data), data)
  # Extract response and predictors from formula
  mf <- model.frame(formula,data)
  assign(as.character(formula[[2]]), model.response(mf))
  environment(formula) <- environment()
  requireNamespace("vegan", quietly = TRUE)
  vegan::adonis2(formula, data)
}
