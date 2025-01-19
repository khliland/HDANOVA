#' Linear Mixed Model PCA
#'
#' @param formula Model formula accepting a single response (block) and predictors. See Details for more information.
#' @param data The data set to analyse.
#' @param pca.in Compress response before ASCA (number of components), default = 5.
#' @param aug_error Error term of model ("denominator", "residual", numeric alpha-value).
#' The latter implies the first with a scaling factor.
#' @param use_ED Use Effective Dimensions instead of degrees of freedom when scaling.
#' @param REML Use restricted maximum likelihood estimation.
#' Alternatives: TRUE (default), FALSE (ML), NULL (least squares).
#' @param contrasts Effect coding: "sum" (default = sum-coding), "weighted", "reference", "treatment".
#' @param ... Additional arguments to \code{\link{asca_fit}}.
#'
#' @description
#' This function mimics parts of the LiMM-PCA framework, combining ASCA+ and
#' linear mixed models to analyse high-dimensional designed data. The default
#' is to use REML estimation and scaling of the backprojected errors. See examples
#' for alternatives.
#'
#' @references
#' * Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and linear mixed models to analyse high-dimensional designed data. Journal of Chemometrics, 34(6), e3232.
#'
#' @return An object of class \code{limmpca}, inheriting from the general \code{asca} class.
#' @export
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{asca_fit}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#'
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Default LiMM-PCA model with two factors and interaction, 5 PCA components
#' mod <- limmpca(assessment ~ candy*r(assessor), data=candies)
#' summary(mod)
#' scoreplot(mod, factor = "candy")
#'
#' # LiMM-PCA with least squares estimation and 8 PCA components
#' modLS <- limmpca(assessment ~ candy*r(assessor), data=candies, REML=NULL, pca.in=8)
#' summary(modLS)
#' scoreplot(modLS, factor = "candy")
#'
#' # Load Caldana data
#' data(caldana)
#'
#' # Combining effects in LiMM-PCA (assuming light is a random factor)
#' mod.comb <- limmpca(compounds ~ time + comb(r(light) + r(time:light)), data=caldana, pca.in=8)
#' summary(mod.comb)
limmpca <- function(formula, data, pca.in = 5, aug_error = 0.05,
                    use_ED = FALSE, REML = TRUE, contrasts = "contr.sum", ...){
  # formula, data, subset, weights, na.action, family, permute=FALSE,
  # unrestricted = FALSE,
  # add_error = FALSE, # TRUE => APCA/LiMM-PCA
  # aug_error = "denominator", # "residual" => Mixed, alpha-value => LiMM-PCA
  # pca.in = FALSE, # n>1 => LiMM-PCA and PC-ANOVA
  # coding = c("sum","weighted","reference","treatment"),
  # SStype = "II",
  # REML = NULL
  object <- asca_fit(formula=formula, data=data, pca.in=pca.in,
                     aug_error=aug_error, use_ED=use_ED, REML=REML,
                     contrasts = contrasts, ...)
  object$call <- match.call()
  class(object) <- c("limmpca", class(object))
  object
}
