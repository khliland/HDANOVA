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
#' # Result plotting for first factor
#' loadingplot(mod, scatter=TRUE, labels="names")
#' scoreplot(mod)
#' # ... and second factor
#' scoreplot(mod, factor = "assessor")
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
pcanova <- function(formula, data, ncomp, subset, weights, na.action, family, SStype = "III"){
  ## Force contrast to sum
  opt <- options(contrasts = c(unordered="contr.sum", ordered="contr.poly"))
  on.exit(options(opt))

  ## Get the data matrices
  Y <- data[[formula[[2]]]]
  N <- nrow(Y)
  p <- ncol(Y)
  Y <- Y - rep(colMeans(Y), each=N) # Centre Y
  ssqY <- sum(Y^2)
  if(missing(ncomp)){
    # Determine automatically by cross-validation
    ncomp <- which.min(.PCAcv(Y))
  } else {
    if(ncomp < 1 && ncomp > 0){
      # Find first component with cummulative explained variance higher than or equal to selected proportion
      pca <- .pca(Y)
      ncomp <- min(which(cumsum(pca$explvar/100) >= ncomp))
    }
  }
  pca <- .pca(Y, ncomp = ncomp)

  # Yudv <- svd(Y)
  #Y <- Yudv$u[,1:pca.in,drop=FALSE] * rep(Yudv$d[1:pca.in], each=N)
  Y <- pca$scores[, 1:ncomp, drop=FALSE]

  # Prepare for storage
  aovs <- list()
  mods <- list()

  mf <- match.call(expand.dots = FALSE)
  fit.type <- "'lm' (Linear Model)"
  if(length(grep('|', formula, fixed=TRUE)) == 0){
    # Fixed effect model
    if(missing(family)){
      # LM
      m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("lm")
      mf[[3]] <- as.name("dat")
      dat <- data
      for(i in 1:ncol(Y)){
        dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
        ano   <- eval(mf, envir = environment())
        mods[[i]] <- ano
        if(i == 1)
          coefs <- matrix(0.0, length(coefficients(ano)), ncol(Y))
        coefs[,i] <- coefficients(ano)
        if(SStype == "I" || SStype == 1)
          aovs[[i]] <- anova(ano)
        if(SStype == "II" || SStype == 2)
          aovs[[i]] <- Anova(ano, type = "II")
        if(SStype == "III" || SStype == 3)
          aovs[[i]] <- Anova(ano, type = "III")
      }
    } else {
      # GLM
      m <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("glm")
      mf[[3]] <- as.name("dat")
      dat <- data
      for(i in 1:ncol(Y)){
        dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
        ano <- eval(mf, envir = environment())
        mods[[i]] <- ano
        if(i == 1)
          coefs <- matrix(0.0, length(coefficients(ano)), ncol(Y))
        coefs[,i] <- coefficients(ano)
        if(SStype == "I" || SStype == 1)
          aovs[[i]] <- anova(ano)
        if(SStype == "II" || SStype == 2)
          aovs[[i]] <- Anova(ano, type = "II")
        if(SStype == "III" || SStype == 3)
          aovs[[i]] <- Anova(ano, type = "III")
      }
      fit.type <- "'glm' (Generalized Linear Model)"
    }
  } else {
    # Mixed model
    if(missing(family)){
      # LM
      m <- match(c("formula", "data", "weights", "subset", "na.action"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("lmer")
      mf[[3]] <- as.name("dat")
      dat <- data
      for(i in 1:ncol(Y)){
        dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
        ano <- eval(mf, envir = environment())
        mods[[i]] <- ano
        if(i == 1)
          coefs <- matrix(0.0, length(colMeans(coefficients(ano)[[1]])), ncol(Y))
        coefs[,i] <- colMeans(coefficients(ano)[[1]])
        if(SStype == "I" || SStype == 1)
          aovs[[i]] <- anova(ano)
        if(SStype == "II" || SStype == 2)
          aovs[[i]] <- Anova(ano, type = "II")
        if(SStype == "III" || SStype == 3)
          aovs[[i]] <- Anova(ano, type = "III")
      }
      fit.type <- "'lmer' (Linear Mixed Model)"
    } else {
      # GLM
      m <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
      mf <- mf[c(1, m)]                # Retain only the named arguments
      mf[[1]] <- as.name("glmer")
      mf[[3]] <- as.name("dat")
      dat <- data
      for(i in 1:ncol(Y)){
        dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
        ano <- eval(mf, envir = environment())
        mods[[i]] <- ano
        if(i == 1) # colMeans assumes only random intercepts, not slopes
          coefs <- matrix(0.0, length(colMeans(coefficients(ano)[[1]])), ncol(Y))
        coefs[,i] <- colMeans(coefficients(ano)[[1]])
        if(SStype == "I" || SStype == 1)
          aovs[[i]] <- anova(ano)
        if(SStype == "II" || SStype == 2)
          aovs[[i]] <- Anova(ano, type = "II")
        if(SStype == "III" || SStype == 3)
          aovs[[i]] <- Anova(ano, type = "III")
      }
      fit.type <- "'glmer' (Generalized Linear Mixed Model)"
    }
  }
  names(mods) <- names(aovs) <- paste0("Comp. ", 1:length(aovs))
  # M      <- model.matrix(ano)
  # effs   <- attr(terms(ano), "term.labels")
  # assign <- attr(M, "assign")
  # modFra <- model.frame(ano)
  #
  # # Exclude numeric effects and their interactions
  # nums   <- names(unlist(lapply(modFra, class)))[which(unlist(lapply(modFra, class)) %in% c("numeric","integer"))]
  # if(length(nums)>0){
  #   exclude  <- match(nums, rownames(attr(terms(ano), "factors")))
  #   approved <- which(colSums(attr(terms(ano), "factors")[exclude,,drop=FALSE])==0)
  # } else {
  #   approved <- 1:max(assign)
  # }
  # if(length(approved)==0)
  #   stop('No factors in model')
  #
  # # Effect loop
  # LS <- effects <- ssq <- list()
  # for(i in 1:length(approved)){
  #   a <- approved[i]
  #   LS[[effs[a]]] <- M[, assign==a, drop=FALSE] %*% coefs[assign==a,]
  #   effects[[effs[a]]] <- modFra[[effs[a]]]
  #
  #   if(i == 1){
  #     residuals <- Y - LS[[effs[i]]]
  #     ssq[[effs[a]]] <- sum(LS[[effs[a]]]^2)
  #   } else {
  #     LSseq <- M[, assign%in%approved[1:i], drop=FALSE] %*% coefs[assign%in%approved[1:i],]
  #     residuals <- Y - LSseq
  #     ssq[[effs[a]]] <- sum(LSseq^2)
  #   }
  # }
  # ssq$res <- ssqY
  # ssq <- unlist(ssq)
  # ssq <- c(ssq[1],diff(ssq))
  # # ssq$res <- sum(residuals^2)
  #
  # # SCAs
  # scores <- loadings <- projected <- singulars <- list()
  # for(i in approved){
  #   maxDir <- min(sum(assign==i), p)
  #   if(pca.in != 0)
  #     maxDir <- min(maxDir, pca.in)
  #   udv <- svd(LS[[effs[i]]])
  #   expli <- (udv$d^2/sum(udv$d^2)*100)[1:maxDir]
  #   scores[[effs[i]]]    <- (udv$u * rep(udv$d, each=N))[,1:maxDir, drop=FALSE]
  #   dimnames(scores[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
  #   loadings[[effs[i]]]  <- udv$v[,1:maxDir, drop=FALSE]
  #   dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
  #   projected[[effs[i]]] <- residuals %*% loadings[[effs[i]]]
  #   dimnames(projected[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
  #   singulars[[effs[i]]] <- udv$d[1:maxDir]
  #   names(singulars[[effs[i]]]) <- paste("Comp", 1:maxDir, sep=" ")
  #   attr(scores[[effs[i]]], 'explvar') <- attr(loadings[[effs[i]]], 'explvar') <- attr(projected[[effs[i]]], 'explvar') <- expli
  #   if(pca.in!=0){ # Transform back if PCA on Y has been performed
  #     loadings[[effs[i]]] <- Yudv$v[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
  #     dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
  #   }
  # }

  obj <- pca
  obj$anovas <- aovs
  obj$coefficients <- coefs
  obj$models <- mods
  obj$call <- match.call()
  obj$fit.type <- fit.type
  #  obj <- list(scores=scores, loadings=loadings, projected=projected, singulars=singulars,
  #              LS=LS, effects=effects, coefficients=coefs, Y=Y, X=M, residuals=residuals,
  #              ssq=ssq, ssqY=ssqY, explvar=ssq/ssqY,
  #              call=match.call(), fit.type=fit.type)
  class(obj) <- c('pcanova', 'pca')
  return(obj)
}
