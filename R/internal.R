#' @importFrom mixlm random.worker
#' @importFrom lme4 lmer
.asca <- function(formula, data, subset, weights, na.action, family, permute,
                  pca.in = FALSE, coding = c("sum","weighted","reference","treatment"),
                  SStype = "III", REML = NULL){

  ## Get the data matrices
  Y <- data[[formula[[2]]]]
  N <- nrow(Y)
  p <- ncol(Y)
  Y <- Y - rep(colMeans(Y), each=N) # Centre Y
  ssqY <- sum(Y^2)
  if(pca.in != 0){
    if(pca.in == 1)
      stop('pca.in = 1 is not supported (single response)')
    Yudv <- svd(Y)
    Y <- Yudv$u[,1:pca.in,drop=FALSE] * rep(Yudv$d[1:pca.in], each=N)
  }
  residuals <- Y

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "weights", "subset", "na.action", "family"), names(mf), 0)
  mf <- mf[c(1, m)]                # Retain only the named arguments

  # Check for random effects using mixlm notation r()
  if( any(grepl("r(",formula,fixed=TRUE)) ){
    rw <- mixlm::random.worker(formula, data, REML)
    mixed <- TRUE
  } else {
    rw <- list(0)
    mixed <- FALSE
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
      fit.func <- "lm"
    } else {
      # GLMM
      fit.type <- "'glmer' (Generalized Linear Mixed Model)"
      fit.func <- "glmer" # TODO: Check mixlm and glmer
    }
  }

  # Simplify SStype
  if(is.character(SStype))
    SStype <- nchar(SStype)

  # Pre-run of model to extract useful information and names
  mf[[1]] <- as.name(fit.func)
  mf[[3]] <- as.name("dat")
  dat <- data
  dat[[formula[[2]]]] <- Y[,1,drop=FALSE]
  opt <- options(contrasts = c(unordered="contr.sum", ordered="contr.poly"))
  mod <- eval(mf, envir = environment())
  effs   <- attr(terms(mod), "term.labels")
  M      <- model.matrix(mod)
  assign <- attr(M, "assign")
  options(opt)
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
      contrasts(dat[[a]]) <- contr.treatment(levels(dat_a))
  }

  # Main ANOVA loop over all responses
  mf[[1]] <- as.name(fit.func)
  mf[[3]] <- as.name("dat")
  sel <- c(names(approved), "Residuals")
  ssq <- numeric(length(sel))
  names(ssq) <- sel
  for(i in 1:ncol(Y)){
    dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
    mod <- eval(mf, envir = environment())
    if(mixed){
      if(i == 0)
        coefs <- matrix(0.0, length(fixef(mod)), ncol(Y))
      coefs[,i] <- fixef(mod) # Gjelder kun ved lmer/glmer
    } else {
      if(i == 1)
        coefs <- matrix(0.0, length(coefficients(mod)), ncol(Y))
      coefs[,i] <- coefficients(mod)
    }
    if(SStype == 1)
      ano <- anova(mod)
    if(SStype == 2)
      ano <- car::Anova(mod, type="II")
    if(SStype == 3)
      ano <- car::Anova(mod, type="III")
    ssq <- ssq + ano[sel,"Sum Sq"]
  }

  M      <- model.matrix(mod)
  effs   <- attr(terms(mod), "term.labels")
  assign <- attr(M, "assign")
  modFra <- model.frame(mod)

  # Permutation testing
  if(!missing(permute)){
    for(i in 1:length(approved)){
      a <- approved[i]
      approved_a <- setdiff(approved, a)
      # Permutation vector
      perm <- 1:N
      for(perm in 1:permute){
        # Loop over other factors
        for(j in approved_a){
          for(k in levels()){

          }
        }
      }
    }
  }

  # Effect loop
  LS <- effects <- list() #  <- ssq
  for(i in 1:length(approved)){
    a <- approved[i]
    LS[[effs[a]]] <- M[, assign==a, drop=FALSE] %*% coefs[assign==a,]
    effects[[effs[a]]] <- modFra[[effs[a]]]

    if(i == 1){
      residuals <- Y - LS[[effs[i]]]
      ssq[[effs[a]]] <- sum(LS[[effs[a]]]^2)
    } else {
      LSseq <- M[, assign%in%approved[1:i], drop=FALSE] %*% coefs[assign%in%approved[1:i],]
      residuals <- Y - LSseq
      #ssq[[effs[a]]] <- sum(LSseq^2)
    }
  }
  #ssq$res <- ssqY
  #ssq <- unlist(ssq)
  #ssq <- c(ssq[1],diff(ssq))

  list(coefs, modFra)
}

dummyvar <- function(x){
  dv <- model.matrix(~x-1, data.frame(x=x, check.names=FALSE))
  colnames(dv) <- levels(x)
  dv
}

SS <- function(x) sum(x^2)
