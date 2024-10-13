#' @name asca
#' @aliases asca
#' @title Analysis of Variance Simultaneous Component Analysis - ASCA
#'
#' @param formula Model formula accepting a single response (block) and predictors. See Details for more information.
#' @param data The data set to analyse.
# #' @param subset Expression for subsetting the data before modelling.
# #' @param weights Optional object weights.
# #' @param na.action How to handle NAs (no action implemented).
# #' @param family Error distributions and link function for Generalized Linear Models.
# #' @param permute Perform approximate permutation testing, default = FALSE (numeric or TRUE = 1000).
# #' @param pca.in Compress response before ASCA (number of components).
#' @param ... Additional arguments to \code{\link{asca_fit}}.
#'
#' @return An \code{asca} object containing loadings, scores, explained variances, etc. The object has
#' associated plotting (\code{\link{asca_plots}}) and result (\code{\link{asca_results}}) functions.
#'
#' @description This is a quite general and flexible implementation of ASCA.
#'
#' @details ASCA is a method which decomposes a multivariate response according to one or more design
#' variables. ANOVA is used to split variation into contributions from factors, and PCA is performed
#' on the corresponding least squares estimates, i.e., \code{Y = X1 B1 + X2 B2 + ... + E = T1 P1' + T2 P2' + ... + E}.
#' This version of ASCA encompasses variants of LiMM-PCA, generalized ASCA and covariates ASCA. It includes
#' confidence ellipsoids for the balanced crossed-effect ASCA.
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
#' @importFrom lme4 lmer
#' @importFrom car ellipse dataEllipse
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{asca_fit}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#' @examples
#' # Load candies data
#' data(candies)
#'
#' # Basic ASCA model with two factors
#' mod <- asca(assessment ~ candy + assessor, data=candies)
#' print(mod)
#'
#' # ASCA model with interaction
#' mod <- asca(assessment ~ candy * assessor, data=candies)
#' print(mod)
#'
#' # Result plotting for first factor
#' loadingplot(mod, scatter=TRUE, labels="names")
#' scoreplot(mod)
#' # No backprojection
#' scoreplot(mod, projections=FALSE)
#' # Spider plot
#' scoreplot(mod, spider=TRUE, projections=FALSE)
#'
#' # ASCA model with compressed response using 5 principal components
#' mod.pca <- asca(assessment ~ candy + assessor, data=candies, pca.in=5)
#'
#' # Mixed Model ASCA, random assessor
#' mod.mix <- asca(assessment ~ candy + r(assessor), data=candies)
#' scoreplot(mod.mix)
#'
#' # Load Caldana data
#' data(caldana)
#'
#' # Combining effects in ASCA
#' mod.comb <- asca(compounds ~ time + comb(light + time:light), data=caldana)
#' summary(mod.comb)
#' timeplot(mod.comb, factor="light", time="time", comb=2)
#'
#' # Permutation testing
#' mod.perm <- asca(assessment ~ candy * assessor, data=candies, permute=TRUE)
#' summary(mod.perm)
#'
#' @export
asca <- function(formula, data, ...){
  # formula, data, subset, weights, na.action, family, permute=FALSE,
  # unrestricted = FALSE,
  # add_error = FALSE, # TRUE => APCA/LiMM-PCA
  # aug_error = "denominator", # "residual" => Mixed, alpha-value => LiMM-PCA
  # pca.in = FALSE, # n>1 => LiMM-PCA and PC-ANOVA
  # coding = c("sum","weighted","reference","treatment"),
  # SStype = "II",
  # REML = NULL
  object <- asca_fit(formula=formula, data=data, ...)
  object$call <- match.call()
  object
}

# asca <- function(formula, data, subset, weights, na.action, family, pca.in = FALSE){
#   ## Force contrast to sum
#   opt <- options(contrasts = c(unordered="contr.sum", ordered="contr.poly"))
#   on.exit(options(opt))
#
#   ## Get the data matrices
#   Y <- data[[formula[[2]]]]
#   N <- nrow(Y)
#   p <- ncol(Y)
#   Y <- Y - rep(colMeans(Y), each=N) # Centre Y
#   ssqY <- sum(Y^2)
#   if(pca.in != 0){
#     if(pca.in == 1)
#       stop('pca.in = 1 is not supported (single response)')
#     Yudv <- svd(Y)
#     Y <- Yudv$u[,1:pca.in,drop=FALSE] * rep(Yudv$d[1:pca.in], each=N)
#   }
#   residuals <- Y
#
#   mf <- match.call(expand.dots = FALSE)
#   fit.type <- "'lm' (Linear Model)"
#   if(length(grep('|', formula, fixed=TRUE)) == 0){
#     # Fixed effect model
#     if(missing(family)){
#       # LM
#       m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
#       mf <- mf[c(1, m)]                # Retain only the named arguments
#       mf[[1]] <- as.name("lm")
#       mf[[3]] <- as.name("dat")
#       dat <- data
#       dat[[formula[[2]]]] <- Y
#       ano   <- eval(mf, envir = environment())
#       coefs <- as.matrix(coefficients(ano))
#     } else {
#       # GLM
#       m <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
#       mf <- mf[c(1, m)]                # Retain only the named arguments
#       mf[[1]] <- as.name("glm")
#       mf[[3]] <- as.name("dat")
#       dat <- data
#       for(i in 1:ncol(Y)){
#         dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
#         ano <- eval(mf, envir = environment())
#         if(i == 1)
#           coefs <- matrix(0.0, length(coefficients(ano)), ncol(Y))
#         coefs[,i] <- coefficients(ano)
#       }
#       fit.type <- "'glm' (Generalized Linear Model)"
#     }
#   } else {
#     # Mixed model
#     if(missing(family)){
#       # LM
#       m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
#       mf <- mf[c(1, m)]                # Retain only the named arguments
#       mf[[1]] <- as.name("lmer")
#       mf[[3]] <- as.name("dat")
#       dat <- data
#       for(i in 1:ncol(Y)){
#         dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
#         ano <- eval(mf, envir = environment())
#         if(i == 1)
#           coefs <- matrix(0.0, length(colMeans(coefficients(ano)[[1]])), ncol(Y))
#         coefs[,i] <- colMeans(coefficients(ano)[[1]])
#       }
#       fit.type <- "'lmer' (Linear Mixed Model)"
#     } else {
#       # GLM
#       m <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
#       mf <- mf[c(1, m)]                # Retain only the named arguments
#       mf[[1]] <- as.name("glmer")
#       mf[[3]] <- as.name("dat")
#       dat <- data
#       for(i in 1:ncol(Y)){
#         dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
#         ano <- eval(mf, envir = environment())
#         if(i == 1) # colMeans assumes only random intercepts, not slopes
#           coefs <- matrix(0.0, length(colMeans(coefficients(ano)[[1]])), ncol(Y))
#         coefs[,i] <- colMeans(coefficients(ano)[[1]])
#       }
#       fit.type <- "'glmer' (Generalized Linear Mixed Model)"
#     }
#   }
#   M      <- model.matrix(ano)
#   effs   <- attr(terms(ano), "term.labels")
#   assign <- attr(M, "assign")
#   modFra <- model.frame(ano)
#
#   # Exclude numeric effects and their interactions
#   nums   <- names(unlist(lapply(modFra, class)))[which(unlist(lapply(modFra, class)) %in% c("numeric","integer"))]
#   if(length(nums)>0){
#     exclude  <- match(nums, rownames(attr(terms(ano), "factors")))
#     approved <- which(colSums(attr(terms(ano), "factors")[exclude,,drop=FALSE])==0)
#   } else {
#     approved <- 1:max(assign)
#   }
#   if(length(approved)==0)
#     stop('No factors in model')
#
#   # Effect loop
#   LS <- effects <- ssq <- list()
#   for(i in 1:length(approved)){
#     a <- approved[i]
#     LS[[effs[a]]] <- M[, assign==a, drop=FALSE] %*% coefs[assign==a,]
#     effects[[effs[a]]] <- modFra[[effs[a]]]
#
#     if(i == 1){
#       residuals <- Y - LS[[effs[i]]]
#       ssq[[effs[a]]] <- sum(LS[[effs[a]]]^2)
#     } else {
#       LSseq <- M[, assign%in%approved[1:i], drop=FALSE] %*% coefs[assign%in%approved[1:i],]
#       residuals <- Y - LSseq
#       ssq[[effs[a]]] <- sum(LSseq^2)
#     }
#   }
#   ssq$res <- ssqY
#   ssq <- unlist(ssq)
#   ssq <- c(ssq[1],diff(ssq))
#   # ssq$res <- sum(residuals^2)
#
#   # SCAs
#   scores <- loadings <- projected <- singulars <- list()
#   for(i in approved){
#     maxDir <- min(sum(assign==i), p)
#     if(pca.in != 0)
#       maxDir <- min(maxDir, pca.in)
#     udv <- svd(LS[[effs[i]]])
#     expli <- (udv$d^2/sum(udv$d^2)*100)[1:maxDir]
#     scores[[effs[i]]]    <- (udv$u * rep(udv$d, each=N))[,1:maxDir, drop=FALSE]
#     dimnames(scores[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
#     loadings[[effs[i]]]  <- udv$v[,1:maxDir, drop=FALSE]
#     dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
#     projected[[effs[i]]] <- residuals %*% loadings[[effs[i]]]
#     dimnames(projected[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
#     singulars[[effs[i]]] <- udv$d[1:maxDir]
#     names(singulars[[effs[i]]]) <- paste("Comp", 1:maxDir, sep=" ")
#     attr(scores[[effs[i]]], 'explvar') <- attr(loadings[[effs[i]]], 'explvar') <- attr(projected[[effs[i]]], 'explvar') <- expli
#     if(pca.in!=0){ # Transform back if PCA on Y has been performed
#       loadings[[effs[i]]] <- Yudv$v[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
#       dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
#     }
#   }
#
#   obj <- list(scores=scores, loadings=loadings, projected=projected, singulars=singulars,
#               LS=LS, effects=effects, coefficients=coefs, Y=Y, X=M, residuals=residuals,
#               ssq=ssq, ssqY=ssqY, explvar=ssq/ssqY,
#               call=match.call(), fit.type=fit.type)
#   if(pca.in!=0){
#     obj$Ypca <- list(svd=Yudv, ncomp=pca.in)
#   }
#   class(obj) <- c('asca', 'list')
#   return(obj)
#   # # Experimental features
#   # # Generalised ASCA, here with a mock Gaussian distribution
#   # mod.glm <- asca(y~x+z, data=dataset, family="gaussian")
#   #
#   # # Generalised Mixed Model ASCA
#   # mod <- asca(y~x+(1|z), data=dataset, family="gaussian")
# }
