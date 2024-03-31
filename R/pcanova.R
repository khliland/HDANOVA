# TODO:
# - Add scoreplot and loadingplot functions
#' @name pcanova
#' @aliases pcanova
#' @title Principal Components Analysis of Variance Simultaneous Component Analysis - PC-ANOVA
#'
#' @param formula Model formula accepting a single response (block) and predictor names separated by + signs.
#' @param data The data set to analyse.
#' @param ncomp The number of components to retain, proportion of variation or default = minimum cross-validation error.
#' @param subset Expression for subsetting the data before modelling.
#' @param weights Optional object weights.
#' @param subset Subset of objects
#' @param na.action How to handle NAs (no action implemented).
#' @param family Error distributions and link function for Generalized Linear Models.
#' @param SStype Sums-of-squares type for Analysis of Variance (I/II/III), defaults to "III".
#'
#' @return TODO: Update. A \code{pcanova} object containing loadings, scores, explained variances, etc. The object has
#' associated plotting (\code{\link{pcanova_plots}}) and result (\code{\link{pcanova_results}}) functions.
#'
#' @description This is a quite general and flexible implementation of PC-ANOVA.
#'
#' @details TODO: Update details. PC-ANOVA is a method which decomposes a multivariate response according to one or more design
#' variables. ANOVA is used to split variation into contributions from factors, and PCA is performed
#' on the corresponding least squares estimates, i.e., \code{Y = X1 B1 + X2 B2 + ... + E = T1 P1' + T2 P2' + ... + E}.
#' This version of PC-ANOVA encompasses variants of LiMM-PCA, generalized PC-ANOVA and covariates PC-ANOVA. It includes
#' confidence ellipsoids for the balanced fixed effect PC-ANOVA.
#'
#' @references Luciano G, Næs T. Interpreting sensory data by combining principal
#' component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167-175.
#'
#' @importFrom lme4 lmer
#' @importFrom car ellipse dataEllipse Anova
#' @importFrom progress progress_bar
#' @importFrom RSpectra svds
#' @seealso TODO
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic PC-ANOVA model with two factors, cross-validated opt. of #components
#' mod <- pcanova(assessment ~ candy + assessor, data = candies)
#' print(mod)
#'
#' # PC-ANOVA model with interaction, minimum 90% explained variance
#' mod <- pcanova(assessment ~ candy * assessor, data = candies, ncomp = 0.9)
#' print(mod)
#' summary(mod)
#'
#' # Tukey group letters for 'candy' per component
#' if("mixlm" %in% rownames(installed.packages())){
#'   lapply(mod$models, function(x)
#'          mixlm::cld(mixlm::simple.glht(x,
#'                                        effect = "candy")))
#' }
#'
#' # Result plotting
#' loadingplot(mod, scatter=TRUE, labels="names")
#' scoreplot(mod)
#'
#' # Mixed Model PC-ANOVA, random assessor
#' mod.mix <- pcanova(assessment ~ candy + (1|assessor), data=candies, ncomp = 0.9)
#' scoreplot(mod.mix)
#' # Fixed effects
#' summary(mod.mix)
#' # Random effects
#' lapply(mod.mix$models, lme4::ranef)
#'
#' @export
pcanova <- function(formula, data, ncomp = 0.9, ...){
  # Run ASCA
  object <- asca_fit(formula, data, pca.in = ncomp, ...)
  # Extract relevant parts
  pc <- object$Ypca$pca
  pc$anovas <- object$anovas
  pc$fit.type <- object$fit.type
  pc$Y <- object$Y
  pc$X <- object$X
  pc$models <- object$models
  # Remove redundant parts
#  object[c("scores","loadings","projected","singulars","LS","effects","coefficients","residuals","ssq","ssqY","explvar")] <- NULL
  # Explained variance
#  object$explvar <- object$Ypca$svd$d^2/sum(object$Ypca$svd$d^2)
  names(pc$anovas) <- paste0("Comp. ", 1:length(pc$anovas))
  pc$call <- match.call()
  class(pc) <- c('pcanova', 'mvr')
  pc
}

#' @export
#' @rdname pcanova
print.pcanova <- function(x, ...){
  cat("PC-ANOVA - Principal Components Analysis of Variance\n")
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x$anovas)
}

#' @export
#' @rdname pcanova
summary.pcanova <- function(object, ...){
  cat("PC-ANOVA - Principal Components Analysis of Variance\n")
  cat("\nCall:\n", deparse(object$call), "\n", sep = "")
  print(object$anovas)
}


# # Experimental features
# # Generalised PC-ANOVA, here with a mock Gaussian distribution
# mod.glm <- pcanova(y~x+z, data=dataset, family="gaussian")
#
# # Generalised Mixed Model PC-ANOVA
# mod <- pcanova(y~x+(1|z), data=dataset, family="gaussian")
