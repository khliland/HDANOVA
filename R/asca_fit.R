# Weighted coding is only available for lm, not lme4??
# @importFrom mixlm random.worker lm
#' @importFrom lme4 lmer glmer
#' @importFrom progress progress_bar
#' @export
asca_fit <- function(formula, data, subset, weights, na.action, family, permute,
                     add_error = FALSE, # TRUE => APCA/LiMM-PCA
                     aug_error = "residual", # "denominator" => Mixed, alpha-value => LiMM-PCA
                     pca.in = FALSE, # n>1 => LiMM-PCA and PC-ANOVA
                     coding = c("sum","weighted","reference","treatment"),
                     SStype = "III", REML = NULL){

  # Simplify SStype
  if(is.character(SStype))
    SStype <- nchar(SStype)

  ## Get the data matrices
  Yorig <- Y <- data[[formula[[2]]]]
  N <- nrow(Y)
  p <- ncol(Y)
  Y <- Y - rep(colMeans(Y), each=N) # Centre Y
  ssqY <- sum(Y^2)
  if(pca.in != 0){ # Pre-decomposition, e.g., LiMM-PCA, PC-ANOVA
    if(is.numeric(pca.in) && pca.in == 1)
      stop('pca.in = 1 is not supported (single response)')
    # Automatic determination of dimension
    if(is.logical(pca.in) && pca.in)
      pca.in <- which.min(.PCAcv(Y))
    # Dimensions according to explained variance
    if(pca.in < 1){
      pca <- .pca(Y)
      pca.in <- min(which(cumsum(pca$explvar/100) >= pca.in))
    }
    # PCA scores
    pca <- .pca(Y, ncomp = pca.in)
    Y <- pca$scores
    #    Yudv <- svd(Y)
    #    Y <- Yudv$u[,1:pca.in,drop=FALSE] * rep(Yudv$d[1:pca.in], each=N)
  }
  residuals <- Y

  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
  mf <- mf[c(1, m)]                # Retain only the named arguments
  lme4 <- FALSE

  # Check for random effects using lmer notation (lme4) |
  if( any(grepl("|",formula,fixed=TRUE)) ){
    mixed <- TRUE
    rw <- list(rformula = formula)
    lme4 <- TRUE
  } else {
    # Check for random effects using mixlm notation r()
    if( any(grepl("r(",formula,fixed=TRUE)) ){
      rw <- mixlm::random.worker(formula, data, REML)
      #mf$formula <- rw$formula
      mixed <- TRUE
    } else {
      rw <- list(0)
      mixed <- FALSE
    }
  }
  if(!mixed){
    # Fixed effect model
    if(missing(family)){
      # LM
      fit.type <- "'lm' (Linear Model)"
      fit.func <- "lm"
    } else {
      # GLM
      fit.type <- "'glm' (Generalized Linear Model)"
      fit.func <- "glm" # TODO: Check mixlm and glm
    }
  } else {
    # Mixed model
    if(missing(family)){
      # LMM
      fit.type <- "'lmm' (Linear Mixed Model)"
      if(lme4){
        fit.func <- "lmer"
      } else {
        fit.func <- "lm"
      }
    } else {
      # GLMM
      fit.type <- "'glmer' (Generalized Linear Mixed Model)"
      fit.func <- "glmer" # TODO: Check mixlm and glmer
    }
  }

  # Pre-run of model to extract useful information and names
  mfPre <- mf
  mfPre[[1]] <- as.name(fit.func)
  mfPre[[3]] <- as.name("dat")
  dat <- data
  dat[[formula[[2]]]] <- Y[,1,drop=FALSE]
  #opt <- options(contrasts = c(unordered="contr.sum", ordered="contr.poly"))
  mod <- eval(mfPre, envir = environment())
  effs   <- attr(terms(mod), "term.labels")
  M      <- model.matrix(mod)
  assign <- attr(M, "assign")
  #options(opt)
  modFra <- model.frame(mod)

  # Exclude numeric effects and their interactions
  nums <- names(unlist(lapply(modFra, class)))[which(unlist(lapply(modFra, class)) %in% c("numeric","integer"))]
  if(length(nums)>0){
    exclude  <- match(nums, rownames(attr(terms(modFra), "factors")))
    approved <- which(colSums(attr(terms(modFra), "factors")[exclude,,drop=FALSE])==0)
  } else {
    approved <- 1:max(assign)
  }
  if(length(approved)==0)
    stop('No factors in model')
  approvedMain <- which(colSums(attr(terms(mod), "factors"))[approved]==1)
  names(approved) <- effs[approved]

  # Apply coding to all included factors
  if(length(coding)>1)
    coding <- coding[1]
  if(!coding %in% c("sum","weighted","reference","treatment"))
    stop("Invalid coding")
  for(i in 1:length(approvedMain)){
    a <- approvedMain[i]
    dat_a <- dat[[effs[a]]]
    if(!missing(subset))
      dat_a <- subset(dat_a, subset)
    if(coding == "sum")
      contrasts(dat[[effs[a]]]) <- contr.sum(levels(dat_a))
    if(coding == "weighted"){
      contrasts(dat[[effs[a]]]) <- contr.weighted(dat_a)
    }
    if(coding == "reference" || coding == "treatment")
      contrasts(dat[[effs[a]]]) <- contr.treatment(levels(dat_a))
  }
  if(fit.func == "lm"){
    if(coding == "sum")
      mf$contrasts <- mfPre$contrasts <- "contr.sum"
    if(coding == "treatment" || coding == "reference")
      mf$contrasts <- mfPre$contrasts <- "contr.treatment"
    if(coding == "weighted")
      mf$contrasts <- mfPre$contrasts <- "contr.weighted"
  }

  # Main ANOVA loop over all responses
  mf[[1]] <- as.name(fit.func)
  mfPre[[3]] <- as.name("dat")
  sel <- c(names(approved), "Residuals")
  ssq <- numeric(length(sel))
  names(ssq) <- sel
  mod <- eval(mfPre, envir = environment())#, envir = parent.frame())

  #mod <- eval(mf, envir = environment())#, envir = parent.frame())
  anovas <- models <- list()
  if(mixed)
    formulaStr <- as.character(rw$rformula)
  else
    formulaStr <- as.character(formula)
  for(i in 1:ncol(Y)){
    if(mixed){
      dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
#      modi <- update(mod)
      modi <- eval(mod, envir = environment())
      #      modi <- update(mod, formula(paste0("`",formulaStr[2], "[,",i,"]`", formulaStr[1], formulaStr[3])))
      if(lme4){
        if(i == 1)
          coefs <- matrix(0.0, length(lme4::fixef(modi)), ncol(Y))
        coefs[,i] <- lme4::fixef(modi)
      } else {
        if(i == 1)
          coefs <- matrix(0.0, length(coefficients(modi)), ncol(Y))
        coefs[,i] <- coefficients(modi)
        #  coefs <- matrix(0.0, length(fixef(modi)), ncol(Y))
        #coefs[,i] <- fixef(modi) # Gjelder kun ved lmer/glmer
      }
    } else {
      dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
#      modi <- update(mod)
      modi <- eval(mod, envir = environment())
      #      modi <- update(mod, formula(paste0(formulaStr[2], "[,",i,"]", formulaStr[1], formulaStr[3])))
      if(i == 1)
        coefs <- matrix(0.0, length(coefficients(modi)), ncol(Y))
      coefs[,i] <- coefficients(modi)
    }

    if(SStype == 1)
      ano <- anova(modi)
    if(SStype == 2)
      ano <- Anova(modi, type="II")
    if(SStype == 3)
      ano <- Anova(modi, type="III")
    if(mixed)
      ssq <- ssq + ano$anova[sel,"Sum Sq"]
    else
      ssq <- ssq + ano[sel,"Sum Sq"]
    models[[i]] <- modi
    anovas[[i]] <- ano
  }

  M      <- model.matrix(mod)
  effs   <- attr(terms(mod), "term.labels")
  assign <- attr(M, "assign")
  modFra <- HDANOVA::extended.model.frame(model.frame(mod), data)

  # Effect loop
  LS <- effects <- list() #  <- ssq
  for(i in 1:length(approved)){
    a <- approved[i]
    LS[[effs[a]]] <- M[, assign==a, drop=FALSE] %*% coefs[assign==a,]
    effects[[effs[a]]] <- modFra[[effs[a]]]
  }
  # Residuals
  residuals <- Y - M%*%coefs

  # Augment error term to LS for permutation testing, LiMMPCA and similar
  LS_aug <- LS
  if(aug_error == "residuals" || !mixed){ # Fixed effect models and forced "residuals"
    # Add residuals to all LSs (augmented for LiMM-PCA and similar)
    for(i in 1:length(approved)){
      a <- approved[i]
      LS_aug[[effs[a]]] <- LS[[effs[a]]] + residuals
    }
  } else {
    # Augment errors according to Mixed Model ANOVA
    for(i in 1:length(approved)){
      a <- approved[i]
      C <- 1
      browser()
      if(is.numeric(add_error) || add_error){
        # LiMM-PCA -> sqrt(dfNum / dfDenom * F(dfNum, dfDenom, 1-alpha))
        #             sqrt(dfNum / dfDenom * pf(1-alpha, dfNum, dfDenom))
        ets <- ano$err.terms
        names(ets) <- rownames(ano$anova)
        dfDenom <- ano$denom.df
        dfNum <- ano$anova[["Df"]]
        names(dfNum) <- rownames(ano$anova)

        # Effective Dimensions som andel forklart varians forklart av
        # korrekt del av LS-matriser

        alpha <- ifelse(is.numeric(add_error), add_error, 0.05)
        browser()
        C <- sqrt()
      }
    }
  }

  if(add_error){
    # Add residuals to all LSs
    for(i in 1:length(approved)){
      a <- approved[i]
      LS[[effs[a]]] <- LS_aug[[effs[a]]]
    }
  }

  # Permutation testing
  if(!missing(permute)){
    # Default to 1000 permutations
    if(is.logical(permute)){
      permute <- 1000
      cat("Defaulting to 1000 permutations\n")
    }
    ssqa <- pvals <- numeric(length(effs))
    ssqaperm <- lapply(1:length(effs),function(i)"NA")
    names(ssqa) <- names(ssqaperm) <- names(pvals) <- effs
    for(i in 1:length(approved)){
      a <- approved[i]
      perms <- numeric(permute)
      # Subset of design matrix for effect a
      D <- M[, assign==a, drop=FALSE]
      # Base ssq
      ssqa[effs[a]] <- norm(D %*% pracma::pinv(D) %*% LS_aug[[effs[a]]], "F")^2
      pb <- progress_bar$new(total = permute, format = paste0("  Permuting ", effs[a], " (", i,"/",length(approved),") [:bar] :percent (:eta)"))
      for(perm in 1:permute){
        perms[perm] <- norm(D %*% pracma::pinv(D) %*% LS_aug[[effs[a]]][sample(N),], "F")^2
        pb$tick()
      }
      ssqaperm[[effs[a]]] <- perms
      pvals[effs[a]] <- sum(perms > ssqa[effs[a]])/permute
    }
  }

  # SCAs
  scores <- loadings <- projected <- singulars <- list()
  for(i in approved){
    maxDir <- min(sum(assign==i), p)
    if(pca.in != 0)
      maxDir <- min(maxDir, pca.in)
    udv <- svd(LS[[effs[i]]])
    expli <- (udv$d^2/sum(udv$d^2)*100)[1:maxDir]
    scores[[effs[i]]]    <- (udv$u * rep(udv$d, each=N))[,1:maxDir, drop=FALSE]
    dimnames(scores[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
    loadings[[effs[i]]]  <- udv$v[,1:maxDir, drop=FALSE]
    dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
    projected[[effs[i]]] <- residuals %*% loadings[[effs[i]]]
    dimnames(projected[[effs[i]]]) <- list(rownames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
    singulars[[effs[i]]] <- udv$d[1:maxDir]
    names(singulars[[effs[i]]]) <- paste("Comp", 1:maxDir, sep=" ")
    attr(scores[[effs[i]]], 'explvar') <- attr(loadings[[effs[i]]], 'explvar') <- attr(projected[[effs[i]]], 'explvar') <- expli
    if(pca.in!=0){ # Transform back if PCA on Y has been performed
      loadings[[effs[i]]] <- pca$loadings[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
      #loadings[[effs[i]]] <- Yudv$v[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
      dimnames(loadings[[effs[i]]]) <- list(colnames(LS[[effs[i]]]), paste("Comp", 1:maxDir, sep=" "))
    }
  }
  obj <- list(scores=scores, loadings=loadings, projected=projected, singulars=singulars,
              LS=LS, effects=effects, coefficients=coefs, Y=Yorig, X=M, residuals=residuals,
              ssq=ssq, ssqY=ssqY, explvar=ssq/ssqY, models=models, anovas=anovas,
              call=match.call(), fit.type=fit.type)
  if(pca.in!=0){
    obj$Ypca <- list(pca=pca, ncomp=pca.in)
  }
  if(!missing(permute)){
    obj$permute <- list(ssqa=ssqa, ssqaperm=ssqaperm, pvalues=pvals)
  }
  class(obj) <- c('asca', 'list')
  return(obj)
}

dummyvar <- function(x){
  dv <- model.matrix(~x-1, data.frame(x=x, check.names=FALSE))
  colnames(dv) <- levels(x)
  dv
}

SS <- function(x) sum(x^2)


formula_comb <- function(formula, data, REML = NULL){
  formula <- formula(formula)
  terms <- terms(formula)
  effsr <- attr(terms,"term.labels")
  effs  <- attr(terms(cparse(formula)),"term.labels")
  if(length(effs)==0){
    return( list(0) )
  }

  has.intercept <- attr(terms,"intercept")==1
  combs <- sort(unique(c(grep("[:]comb[(]",effsr),   # Match combination in interaction
                         grep("^comb[(]",  effsr),   # Match combination in the beginning
                         grep("[(]comb[(]",effsr))))  # Match combination inside function

  eff.splits <- list()
  for(i in 1:length(effs)){ # Split effect to look for hidden combination interactions
    eff.splits[[i]] <- fparse(formula(paste("1~", effs[i],sep="")))
  }
  eff.lengths <- lapply(eff.splits,length)
  main.effs   <- effs[eff.lengths==1]
  main.combs  <- main.effs[main.effs%in%effs[combs]]
  main.combs.only.inter <- character(0)
  for(i in combs){
    main.combs.only.inter <- c(main.combs.only.inter, setdiff(eff.splits[[i]],main.effs)) # Combination main effects only present in interactions
  }
  inter.combs <- which(unlist(lapply(eff.splits,function(i) any(main.combs%in%i))))
  # Check if any interactions containing combination effects are not labeled as combination
  if(any(is.na(match(inter.combs,combs)))){
    extra.randoms <- inter.combs[which(is.na(match(inter.combs,combs)))]
    warning(paste(paste(effs[extra.randoms],sep="",collapse=", "), " included as combination interaction",ifelse(length(extra.randoms)==1,"","s"),sep=""))
    combs <- cbind(combs,extra.randoms)
    effs  <- effs[!(extra.randoms%in%effs)]
  }
  if(length(combs)==0){
    return( list(0) )
  } else {
    if(is.logical(REML)){
      remleffs     <- c(effs[setdiff(1:length(effs),combs)],paste("(1|",effs[combs],")",sep=""))
      reml.formula <- formula(paste(formula[[2]],"~",paste(remleffs,collapse="+"),ifelse(has.intercept,"","-1"),sep=""))

      return( list(cformula = formula,
                   formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")),
                   combination = effs[combs],
                   main.combs.only.inter = main.combs.only.inter,
                   fixed = effs[setdiff(1:length(effsr),combs)],
                   all = effs,
                   allc = effsr,
                   has.intercept = has.intercept,
                   remleffs = remleffs,
                   reml.formula = reml.formula))
    } else {
      return( list(cformula = formula,
                   formula = formula(paste(formula[[2]],"~",paste(effs,collapse="+"),ifelse(has.intercept,"","-1"),sep="")),
                   combination = effs[combs],
                   main.combs.only.inter = main.combs.only.inter,
                   fixed = effs[setdiff(1:length(effsr),combs)],
                   all = effs,
                   allc = effsr,
                   has.intercept = has.intercept))
    }
  }
}

# Remove comb() from formula (possibly convert to lmer formula)
cparse <- function (f, REML = FALSE) {
  if(!inherits(f,'formula'))
    # if (class(f) != "formula")
    stop("'f' must be a formula")
  right <- attr(terms(f),"term.labels") # Let R split short-hand notations first
  if( length(right) == 0){
    return(f)
  }
  n <- length(right)
  result      <- character(n)
  result.REML <- character(n)

  # Main recursive loop extracting effects without comb()
  for(i in 1:n){
    parsecall <- function(x) {
      if(!inherits(x,'call'))
        # if (class(x) != "call")
        stop("'x' must be a call")
      if (length(x[[1]]) == 1) {
        return(x)
      }
      if (length(x[[1]]) == 2 && x[[1]][[1]] == "comb") {
        return(x[[1]][2])
      }
      for (i in 2:length(x[[1]])) x[[1]][i] <- parsecall(x[[1]][i])
      return(x)
    }
    result[[i]] <- as.character(parsecall(formula(paste("~",paste(right[i],sep="+")))[2]))
  }
  f[3] <- formula(paste("~", paste(result, sep="", collapse="+")))[2]

  # Recursive loop adding (1 | x) notation for REML estimation
  if(REML){
    for(i in 1:n){
      parsecall <- function(x) {
        if(!inherits(x,'call'))
          # if (class(x) != "call")
          stop("'x' must be a call")
        if (length(x[[1]]) == 1) {
          return(FALSE)
        }
        if (length(x[[1]]) == 2 && x[[1]][[1]] == "comb") {
          return(TRUE)
        } else {
          tmp <- logical(length(x[[1]]))
          for (j in 2:length(x[[1]])) tmp[j-1] <- parsecall(x[[1]][j])
          return(any(tmp))
        }
      }
      ran <- parsecall(formula(paste("~",paste(right[i],sep="+")))[2])
      result.REML[i] <- ifelse(ran,as.character(formula(paste("~(1 | ",result[[i]],")",sep=""))[2]), result[[i]])
    }
    f[3] <- formula(paste("~", paste(result.REML, sep="", collapse="+")))[2]
  }
  f
}

# Extract variable names from formula
fparse <- function(f) {
  if(!inherits(f,'formula')) stop("'f' must be a formula")
  # if (class(f) != "formula") stop("'f' must be a formula")
  right <- f[3]
  parsecall <- function(x) {
    if(!inherits(x,'call')) stop("'x' must be a call")
    # if (class(x) != "call") stop("'x' must be a call")
    if (length(x[[1]]) == 1) {
      if(is.numeric(x[[1]])) {
        return()
      } else {
        return(as.character(x[[1]]))
      }
    }
    res <- list()
    for (i in 2:length(x[[1]]))
      res[[i-1]] <- parsecall(x[[1]][i])
    return(unlist(res))
  }
  unique(parsecall(right))
}
