#' @name apls
#' @aliases apls
#' @title Analysis of Variance Partial Least Squares - APLS
#'
#' @param formula Model formula accepting a single response (block) and predictors. See Details for more information.
#' @param data The data set to analyse.
#' @param add_error Add error to LS means (default = TRUE).
#' @param contrasts Effect coding: "sum" (default = sum-coding), "weighted", "reference", "treatment".
#' @param permute Number of permutations to perform (default = 1000).
#' @param perm.type Type of permutation to perform, either "approximate" or "exact" (default = "approximate").
# #' @param subset Expression for subsetting the data before modelling.
# #' @param weights Optional object weights.
# #' @param na.action How to handle NAs (no action implemented).
# #' @param family Error distributions and link function for Generalized Linear Models.
# #' @param permute Perform approximate permutation testing, default = FALSE (numeric or TRUE = 1000).
# #' @param pca.in Compress response before APLS (number of components).
#' @param ... Additional arguments to \code{\link{hdanova}}.
#'
#' @return An \code{apls} object containing loadings, scores, explained variances, etc. The object has
#' associated plotting (\code{\link{asca_plots}}) and result (\code{\link{asca_results}}) functions.
#'
#' @description This is a quite general and flexible implementation of APLS.
#'
#' @details APLS is a method which decomposes a multivariate response according to one or more design
#' variables. ANOVA is used to split variation into contributions from factors, and PLS is performed
#' on the corresponding least squares estimates, i.e., \code{Y = X1 B1 + X2 B2 + ... + E = T1 P1' + T2 P2' + ... + E}.
#' For balanced designs, the PLS components are equivalent to PCA components, i.e., APLS and APCA are equivalent.
#' This version of APLS encompasses variants of LiMM-PLS, generalized APLS and covariates APLS.
#'
#' The formula interface is extended with the function r() to indicate random
#' effects and comb() to indicate effects that should be combined. See Examples
#' for use cases.
#'
#' @references
#' * Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J., and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA): A new tool for analyzing designed metabolomics data. Bioinformatics, 21(13), 3043–3048.
#' * Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence ellipsoids for ASCA models based on multivariate regression theory. Journal of Chemometrics, 32(e2990), 1–13.
#' * Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and linear mixed models to analyse high-dimensional designed data. Journal of Chemometrics, 34(6), e3232.
#'
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{hdanova}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic APLS model with two factors
#' mod <- apls(assessment ~ candy + assessor, data=candies)
#' print(mod)
#'
#' # APLS model with interaction
#' mod <- apls(assessment ~ candy * assessor, data=candies)
#' print(mod)
#'
#' # Result plotting for first factor
#' loadingplot(mod, scatter=TRUE, labels="names")
#' scoreplot(mod)
#' # No backprojection
#' scoreplot(mod, projections=FALSE)
#' # Spider plot
#' scoreplot(mod, spider=TRUE)
#'
#' # APLS model with compressed response using 5 principal components
#' mod.pca <- apls(assessment ~ candy + assessor, data=candies, pca.in=5)
#'
#' # Mixed Model APLS, random assessor
#' mod.mix <- apls(assessment ~ candy + r(assessor), data=candies)
#' scoreplot(mod.mix)
#'
#' # Mixed Model APLS, REML estimation
#' mod.mix <- apls(assessment ~ candy + r(assessor), data=candies, REML=TRUE)
#' scoreplot(mod.mix)
#'
#' # Load Caldana data
#' data(caldana)
#'
#' # Combining effects in APLS
#' mod.comb <- apls(compounds ~ time + comb(light + time:light), data=caldana)
#' summary(mod.comb)
#' timeplot(mod.comb, factor="light", time="time", comb=2)
#'
#' # Permutation testing
#' mod.perm <- apls(assessment ~ candy * assessor, data=candies, permute=TRUE)
#' summary(mod.perm)
#'
#' @export
apls <- function(formula, data, add_error = TRUE, contrasts = "contr.sum",
                 permute = FALSE, perm.type=c("approximate","exact"),...){
  # formula, data, subset, weights, na.action, family, permute=FALSE,
  # unrestricted = FALSE,
  # add_error = FALSE, # TRUE => APCA/LiMM-PCA
  # aug_error = "denominator", # "residual" => Mixed, alpha-value => LiMM-PCA
  # pca.in = FALSE, # n>1 => LiMM-PCA and PC-ANOVA
  # coding = c("sum","weighted","reference","treatment"),
  # SStype = "II",
  # REML = NULL
  object <- hdanova(formula=formula, data=data, add_error = add_error, contrasts = contrasts, ...)
  if(!(is.logical(permute) && !permute)){
    # Default to 1000 permutations
    if(is.logical(permute)){
      permute <- 1000
      if(interactive())
        message("Defaulting to 1000 permutations\n")
    }
    object <- permutation(object, permute=permute, perm.type=perm.type)
  }
  object <- pls(object)
  object$call <- match.call()
  object
}
