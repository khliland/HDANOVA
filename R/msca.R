#' @name msca
#' @aliases msca
#' @title Multilevel Simultaneous Component Analysis - MSCA
#'
#' @param formula Model formula accepting a single response (block) and predictors. See Details for more information.
#' @param data The data set to analyse.
#' @param contrasts Effect coding: "sum" (default = sum-coding), "weighted", "reference", "treatment".
#' @param permute Number of permutations to perform (default = 1000).
#' @param perm.type Type of permutation to perform, either "approximate" or "exact" (default = "approximate").
#' @param ... Additional arguments to \code{\link{hdanova}}.
#'
#' @return An \code{asca} object containing loadings, scores, explained variances, etc. The object has
#' associated plotting (\code{\link{asca_plots}}) and result (\code{\link{asca_results}}) functions.
#'
#' @description This MSCA implementation assumes a single factor to be used as between-individuals factor.
#'
#' @references
#' * Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J., and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA): A new tool for analyzing designed metabolomics data. Bioinformatics, 21(13), 3043–3048.
#' * Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence ellipsoids for ASCA models based on multivariate regression theory. Journal of Chemometrics, 32(e2990), 1–13.
#'
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{hdanova}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic MSCA model with a single factor
#' mod <- msca(assessment ~ candy, data=candies)
#' print(mod)
#' summary(mod)
#'
#' # Result plotting for first factor
#' loadingplot(mod, scatter=TRUE, labels="names")
#' scoreplot(mod)
#'
#' # Within scores
#' scoreplot(mod, factor="within")
#'
#' # Within scores per factor level
#' par.old <- par(mfrow=c(3,2), mar=c(4,4,2,1), mgp=c(2,0.7,0))
#' for(i in 1:length(mod$scores.within))
#'   scoreplot(mod, factor="within", within_level=i,
#'             main=paste0("Level: ", names(mod$scores.within)[i]),
#'             panel.first=abline(v=0,h=0,col="gray",lty=2))
#' par(par.old)
#'
#' # Permutation testing
#' mod.perm <- asca(assessment ~ candy * assessor, data=candies, permute=TRUE)
#' summary(mod.perm)
#'
#' @export
msca <- function(formula, data, contrasts = "contr.sum",
                 permute = FALSE, perm.type=c("approximate","exact"), ...){
  # formula, data, subset, weights, na.action, family, permute=FALSE,
  # unrestricted = FALSE,
  # add_error = FALSE, # TRUE => APCA/LiMM-PCA
  # aug_error = "denominator", # "residual" => Mixed, alpha-value => LiMM-PCA
  # pca.in = FALSE, # n>1 => LiMM-PCA and PC-ANOVA
  # coding = c("sum","weighted","reference","treatment"),
  # SStype = "II",
  # REML = NULL
  object <- hdanova(formula=formula, data=data, contrasts = contrasts, ...)
  if(!(is.logical(permute) && !permute)){
    # Default to 1000 permutations             ------------- Flytt testing til ASCA og venner
    if(is.logical(permute)){
      permute <- 1000
      if(interactive())
        message("Defaulting to 1000 permutations\n")
    }
    object <- permutation(object, permute=permute, perm.type=perm.type)
  }
  object <- sca(object)
  object$call <- match.call()

  # Split within scores per factor level
  between <- object$model.frame[[2]]
  object$scores.within <- list()
  object$explvar.within <- matrix(0, nlevels(between), ncol(object$scores[[2]]))
  for(i in 1:nlevels(between)){
    object$scores.within[[i]] <- object$scores[[2]][between==levels(between)[i],]
    for(j in 1:ncol(object$scores[[2]]))
      object$explvar.within[i,j] <- 100*(1-sum((object$residuals[between==levels(between)[i],,drop=FALSE] - tcrossprod(object$scores[[2]][between==levels(between)[i],j],object$loadings[[2]][,j]))^2)/sum(object$residuals[between==levels(between)[i],,drop=FALSE]^2))
  }
  names(object$scores.within) <- levels(between)
  dimnames(object$explvar.within) <- list(levels(between),colnames(object$scores[[2]]))
  for(i in 1:nlevels(between))
    attr(object$scores.within[[i]], "explvar") <- object$explvar.within[i,]
  class(object) <- c("msca", class(object))
  object
}
