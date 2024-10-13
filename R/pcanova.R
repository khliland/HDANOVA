#' @name pcanova
#' @aliases pcanova
#' @title Principal Components Analysis of Variance Simultaneous Component Analysis - PC-ANOVA
#'
#' @param formula Model formula accepting a single response (block) and predictor names separated by + signs.
#' @param data The data set to analyse.
#' @param ncomp The number of components to retain, proportion of variation or default = minimum cross-validation error.
#' @param ... Additional parameters for the asca_fit function.
#'
#' @return A \code{pcanova} object containing loadings, scores, explained variances, etc. The object has
#' associated plotting (\code{\link{pcanova_plots}}) and result (\code{\link{pcanova_results}}) functions.
#'
#' @description This is a quite general and flexible implementation of PC-ANOVA.
#'
#' @details PC-ANOVA works in the opposite order of ASCA. First the response matrix
#' is decomposed using ANOVA. Then the components are analysed using ANOVA
#' with respect to a design or grouping in the data. The latter can be ordinary
#' fixed effects modelling or mixed models.
#'
#' @references Luciano G, Næs T. Interpreting sensory data by combining principal
#' component analysis and analysis of variance. Food Qual Prefer. 2009;20(3):167-175.
#'
#' @importFrom lme4 lmer
#' @importFrom car ellipse dataEllipse Anova
#' @importFrom progress progress_bar
#' @importFrom RSpectra svds
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{asca_fit}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
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
#' lapply(mod$models, function(x)
#'        mixlm::cld(mixlm::simple.glht(x,
#'                                      effect = "candy")))
#'
#' # Result plotting
#' loadingplot(mod, scatter=TRUE, labels="names")
#' scoreplot(mod)
#'
#' # Mixed Model PC-ANOVA, random assessor
#' mod.mix <- pcanova(assessment ~ candy + r(assessor), data=candies, ncomp = 0.9)
#' scoreplot(mod.mix)
#' # Fixed effects
#' summary(mod.mix)
#'
#' @export
pcanova <- function(formula, data, ncomp = 0.9, ...){
  # Run ASCA
  object <- asca_fit(formula, data, pca.in = ncomp, ...)
  # Extract relevant parts
  pc <- object$Ypca$pca
  pc$Ypca <- object$Ypca
  pc$anovas <- object$anovas
  pc$fit.type <- object$fit.type
  pc$Y <- object$Y
  pc$X <- object$X
  pc$models <- object$models
  pc$effects <- object$effects
  names(pc$anovas) <- paste0("Comp. ", 1:length(pc$anovas))
  pc$call <- match.call()
  class(pc) <- c('pcanova', 'asca', 'list')
  pc
}

# # Experimental features
# # Generalised PC-ANOVA, here with a mock Gaussian distribution
# mod.glm <- pcanova(y~x+z, data=dataset, family="gaussian")
#
# # Generalised Mixed Model PC-ANOVA
# mod <- pcanova(y~x+(1|z), data=dataset, family="gaussian")
