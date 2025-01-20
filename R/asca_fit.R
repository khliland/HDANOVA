# Weighted coding is only available for lm, not lme4??
#' @title ASCA Fitting Workhorse Function
#'
#' @description This function is called by all ASCA related methods in this package. It is documented
#' so that one can have access to a richer set of parameters from the various methods or call this
#' function directly. The latter should be done with care as there are many possibilities and not all
#' have been used in publications or tested thoroughly.
#'
#' @param formula Model formula accepting a single response (block) and predictors. See Details for more information.
#' @param data The data set to analyse.
#' @param subset Expression for subsetting the data before modelling.
#' @param weights Optional object weights.
#' @param na.action How to handle NAs (no action implemented).
#' @param family Error distributions and link function for Generalized Linear Models.
#' @param permute Perform approximate permutation testing, default = FALSE (numeric or TRUE = 1000 permutations).
#' @param perm.type Type of permutation: "approximate" (default) or "exact".
#' @param unrestricted Use unrestricted ANOVA decomposition (default = FALSE).
#' @param add_error Add error to LS means, e.g., for APCA.
#' @param aug_error Augment score matrices in backprojection. Default = "denominator"
#' (of F test), "residual" (force error term), nueric value (alpha-value in LiMM-PCA).
#' @param use_ED Use "effective dimensions" for score rescaling in LiMM-PCA.
#' @param pca.in Compress response before ASCA (number of components).
#' @param contrasts Effect coding: "sum" (default = sum-coding), "weighted", "reference", "treatment".
#' @param coding Defunct. Use 'contrasts' instead.
#' @param equal_baseline Experimental: Set to \code{TRUE} to let interactions, where a main effect is missing,
#' e.g., a nested model, be handled with the same baseline as a cross effect model. If \code{TRUE} the corresponding
#' interactions will be put in quotation marks and included in the \code{model.frame}.
#' @param SStype Type of sum-of-squares: "I" = sequential, "II" (default) = last term, obeying marginality,
#' "III" = last term, not obeying marginality.
#' @param REML Parameter to mixlm: NULL (default) = sum-of-squares, TRUE = REML, FALSE = ML.
#'
#' @return An \code{asca} object containing loadings, scores, explained variances, etc. The object has
#' associated plotting (\code{\link{asca_plots}}) and result (\code{\link{asca_results}}) functions.
#'
#' @importFrom lme4 lmer glmer getME ranef VarCorr fixef
#' @importFrom progress progress_bar
#' @importFrom grDevices adjustcolor palette
#' @importFrom graphics abline axis box hist legend lines points
#' @importFrom stats anova coefficients contr.sum contr.treatment contrasts<- formula getCall model.frame model.matrix model.response qf rnorm sigma terms update
#' @importFrom pracma Rank
#' @export
asca_fit <- function(formula, data, subset, weights, na.action, family,
                     permute=FALSE,
                     perm.type=c("approximate","exact"),
                     unrestricted = FALSE,
                     add_error = FALSE, # TRUE => APCA
                     aug_error = "denominator", # "residual" => Mixed, alpha-value => LiMM-PCA
                     use_ED = FALSE,
                     pca.in = FALSE, # n>1 => LiMM-PCA and PC-ANOVA
                     contrasts = "contr.sum",
                     coding, #c("sum","weighted","reference","treatment"),
                     equal_baseline = FALSE,
                     SStype = "II",
                     REML = NULL){

  # Simplify SStype
  if(is.character(SStype))
    SStype <- nchar(SStype)
  if(!missing(coding))
    stop("Input 'coding' has been exchanged with 'contrasts'. See ?asca_fit.")

  ## Get the data matrices
  #Yorig <- Y <- data[[formula[[2]]]]
  Y <- data[[formula[[2]]]]
  if(inherits(Y, "AsIs"))
    Y <- unclass(Y)
  Yorig <- Y <- as.matrix(Y)
  N <- nrow(Y)
  p <- ncol(Y)
  if(is.null(p))
    stop("Response must be a matrix.")
  #  if(center && (!missing(family) && family!="binomial"))
  #    Y <- Y - rep(colMeans(Y), each=N) # Centre Y by default, but not if family is binomial

  ########################## PCA of response matrix ##########################
  if(pca.in != 0){ # Pre-decomposition, e.g., LiMM-PCA, PC-ANOVA
#    if(is.numeric(pca.in) && pca.in == 1)
#      stop('pca.in = 1 is not supported (single response)')
    # Automatic determination of dimension
    if(is.logical(pca.in) && pca.in)
      pca.in <- which.min(.PCAcv(Y))
    # Dimensions according to explained variance
    if(pca.in < 1){
      pca <- .pca(Y)
      pca.in <- min(which(cumsum(pca$explvar/100) >= pca.in))
    }
    # Limit number of extracted components
    if(pca.in > p){
      warning(paste0("Reducing 'pca.in' from ", pca.in, " to the number of variables (",p,")"))
      pca.in <- p
    }
    # PCA scores
    pca <- .pca(Y, ncomp = pca.in)
    Y <- pca$scores
  }
  ssqY <- sum((Y-rep(colMeans(Y),each=N))^2)
  residuals <- Y

  mf <- match.call(expand.dots = FALSE)
  m  <- match(c("formula", "data", "weights", "subset", "na.action", "family",
                "unrestricted", "REML", "contrasts", "equal_baseline"), names(mf), 0)
  mf <- mf[c(1, m)]                # Retain only the named arguments
  lme4 <- FALSE

  ########################## Decide model type ##########################
  # Check for random effects using lmer notation (lme4) |
  if( any(grepl("|",formula,fixed=TRUE))){
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
  # Add unrestricted to mf for use with mixlm
  if(!lme4 && !is.logical(REML) && is.null(mf$unrestricted))
    mf$unrestricted <- FALSE
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

  ########################## Combined effects ##########################
  # Check for combined factors and remove symbols from formula.
  combined <- FALSE
  if(grepl("comb(", as.character(formula)[3], fixed=TRUE)){
    combined <- list()
    form_no_r <- mixlm::rparse(formula)
    tl <- attr(terms(form_no_r),"term.labels")
    j <- 1
    for(i in 1:length(tl)){
      if(grepl("comb(", tl[i], fixed=TRUE)){
        combined[[j]] <- attr(terms(cparse(formula(paste0(".~", tl[i])))), "term.labels")
        j <- j+1
      }
    }
    formula <- cparse(formula)
  }

  ########################## Dry-run to find properties ##########################
  # Pre-run of model to extract useful information and names
  mfPre <- mf
  mfPre[[1]] <- as.name(fit.func)
  mfPre[[3]] <- as.name("dat")
  dat <- data
  dat[[formula[[2]]]] <- Y[,1,drop=FALSE]
  mod <- eval(mfPre, envir = environment())
  effs   <- attr(terms(mod), "term.labels")
  M      <- model.matrix(mod)
  assign <- attr(M, "assign")
  modFra <- model.frame(mod)

  # Extend effs, assign and M with random effects
  if(lme4 || is.logical(REML)){
    X <- M
    Z <- as.matrix(lme4::getME(mod, "Z"))
    M <- cbind(M, Z)
    u <- lme4::ranef(mod)
    effs <- c(effs, names(u))
    max_assign <- max(assign)
    for(i in 1:length(u)){
      max_assign <- max_assign+1
      assign <- c(assign, rep(max_assign, length(u[[i]][[1]])))
    }
  }

  # Maximum number of dimensions for each effect
  maxDir <- numeric(max(assign))
  for(i in 1:max(assign))
    maxDir[i] <- min(sum(assign==i), p)

  # Alphabetically sorted interactions
  effsAB <- effs
  for(i in 1:length(effs)){
    if(grepl(":", effs[i], fixed=TRUE)){
      effsAB[i] <- paste(sort(strsplit(effs[i],":")[[1]]), collapse=":")
    }
  }

  # Exclude numeric effects and their interactions unless fit.func is lm
  nums <- names(unlist(lapply(modFra, class)))[which(unlist(lapply(modFra, class)) %in% c("numeric","integer"))]
  if(fit.func == "lm")
    nums <- numeric(0)
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
  approvedAB <- approved
  names(approvedAB) <- effsAB[approvedAB]

  ########################## Effect coding ##########################
  # Apply coding to all included factors
  #if(length(coding)>1)
  #  coding <- coding[1]
  #if(!coding %in% c("sum","weighted","reference","treatment"))
  #  stop("Invalid coding")
  # ms <- missing(subset)
  # contrast.list <- lapply(dat, function(dat_a){
  #   if(inherits(dat_a, "factor")){
  #     if(!ms)
  #       dat_a <- subset(dat_a, subset)
  #     if(coding == "sum")
  #       return(contr.sum(levels(dat_a)))
  #     if(coding == "weighted")
  #       return(contr.weighted(dat_a))
  #     if(coding == "reference" || coding == "treatment")
  #       return(contr.treatment(levels(dat_a)))
  #   }
  # })
  # contrast.list <- contrast.list[!sapply(contrast.list, is.null)]

  # for(i in 1:length(approvedMain)){
  #   a <- which(effs==names(approvedMain[i]))
  #   dat_a <- dat[[effs[a]]]
  #   if(!missing(subset))
  #     dat_a <- subset(dat_a, subset)
  #   if(coding == "sum" && is.factor(dat_a))
  #     contrasts(dat[[effs[a]]]) <- contr.sum(levels(dat_a))
  #   if(coding == "weighted" && is.factor(dat_a)){
  #     contrasts(dat[[effs[a]]]) <- contr.weighted(dat_a)
  #   }
  #   if((coding == "reference" || coding == "treatment")  && is.factor(dat_a))
  #     contrasts(dat[[effs[a]]]) <- contr.treatment(levels(dat_a))
  # }
  # if(fit.func == "lm" && !is.logical(REML)){
  #   if(coding == "sum")
  #     mf$contrasts <- mfPre$contrasts <- "contr.sum"
  #   if(coding == "treatment" || coding == "reference")
  #     mf$contrasts <- mfPre$contrasts <- "contr.treatment"
  #   if(coding == "weighted")
  #     mf$contrasts <- mfPre$contrasts <- "contr.weighted"
  # }

  ########################## ANOVA ##########################
  # Main ANOVA loop over all responses
  mf[[1]] <- as.name(fit.func)
  mfPre[[3]] <- as.name("dat")
  sel <- c(effs, "Residuals")
  ssq <- numeric(length(sel))
  names(ssq) <- sel
  mod <- eval(mfPre, envir = environment())
  anovas <- models <- list()
  if(mixed)
    formulaStr <- as.character(rw$rformula)
  else
    formulaStr <- as.character(formula)
  for(i in 1:ncol(Y)){
    if(mixed){
      dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
      modi <- eval(mfPre, envir = environment())
      if(lme4 || is.logical(REML)){
        u <- unlist(lme4::ranef(modi))
        if(i == 1)
          coefs <- matrix(0.0, length(lme4::fixef(modi))+length(u), ncol(Y))
        coefs[,i] <- c(lme4::fixef(modi),u)
      } else {
        if(i == 1)
          coefs <- matrix(0.0, length(coefficients(modi)), ncol(Y))
        coefs[,i] <- coefficients(modi)
      }
    } else {
      dat[[formula[[2]]]] <- Y[,i,drop=FALSE]
      modi <- eval(mfPre, envir = environment())
      if(i == 1){
        coefs <- matrix(0.0, length(coefficients(modi)), ncol(Y))
        rownames(coefs) <- names(coefficients(modi))
      }
      coefs[,i] <- coefficients(modi)
    }

    if(SStype == 1)
      ano <- anova(modi)
    if(SStype == 2){
      if(missing(family))
        ano <- Anova(modi, type="II", singular.ok=TRUE)
      else
        ano <- Anova(modi, type="II", test.statistic = "F", singular.ok=TRUE)
    }
    if(SStype == 3){
      if(missing(family))
        ano <- Anova(modi, type="III", singular.ok=TRUE)
      else
        ano <- Anova(modi, type="III", test.statistic = "F", singular.ok=TRUE)
    }
    if(mixed)
      ssq <- ssq + ano$anova[sel,"Sum Sq"]
    else
      ssq <- ssq + ano[sel,"Sum Sq"]
    models[[i]] <- modi
    anovas[[i]] <- ano
  }
  if(length(ssq) == length(c(effsAB, "Residuals")))
    names(ssq) <- c(effsAB, "Residuals")
  ssq_residual <- ssq[length(ssq)]
  names(models) <- names(anovas) <- colnames(coefs) <- colnames(Y)

  M      <- model.matrix(mod)
  effs   <- attr(terms(mod), "term.labels")
  assign <- attr(M, "assign")
  modFra <- HDANOVA::extended.model.frame(model.frame(mod), data)

  # Extend effs, assign and M with random effects
  if(lme4 || is.logical(REML)){
    M <- cbind(M, as.matrix(getME(mod, "Z")))
    u <- lme4::ranef(mod)
    effs <- c(effs, names(u))
    max_assign <- max(assign)
    for(i in 1:length(u)){
      max_assign <- max_assign+1
      assign <- c(assign, rep(max_assign, length(u[[i]][[1]])))
    }
  }

  # Sort interaction names alphabetically
  for(i in 1:length(colnames(modFra))){
    if(grepl(":", colnames(modFra)[i], fixed=TRUE)){
      colnames(modFra)[i] <- paste(sort(strsplit(colnames(modFra)[i],":")[[1]]), collapse=":")
    }
  }
  # Sort interaction names alphabetically
  for(i in 1:length(effs)){
    if(grepl(":", effs[i], fixed=TRUE)){
      effs[i] <- paste(sort(strsplit(effs[i],":")[[1]]), collapse=":")
    }
  }

  ########################## LS estimates ##########################
  # Effect loop
  LS <- effects <- list()
  for(i in 1:length(approved)){
    a <- approved[i]
    # Exclude non-estimable levels:
    estimable <- !is.na(rowSums(coefs[assign==a,, drop=FALSE]))
    LS[[effs[a]]] <- M[, assign==a, drop=FALSE][,estimable,drop=FALSE] %*% coefs[assign==a,, drop=FALSE][estimable,,drop=FALSE]
    colnames(LS[[effs[a]]]) <- colnames(Y)
    effects[[effs[a]]] <- modFra[[effs[a]]]
  }

  # Residuals
  # Exclude non-estimable levels:
  estimable <- !is.na(rowSums(coefs))
  residuals <- Y - M[,estimable,drop=FALSE] %*% coefs[estimable,,drop=FALSE]

  # If model is of lme4 type, the sums-of-squares are not directly available
  if(length(ssq) == 0){
    for(i in 1:length(approved)){
      a <- approved[i]
      ssq[effs[a]] <- sum(LS[[effs[a]]]^2)
    }
    names(ssq) <- names(effects)
    ssq_residual <- sum(residuals^2)
    ssq[length(ssq)+1] <- ssq_residual
    names(ssq)[length(ssq)] <- "Residuals"
  }

  ########################## Augmented LS/errors ##########################
  # Augment error term to LS for permutation testing, LiMM-PCA and similar
  error <- LS_aug <- LS
  if(aug_error == "residuals" || !mixed){ # Fixed effect models and forced "residuals"
    # Input to asca_fit: aug_error = "residual", # "denominator" => Mixed, alpha-value => LiMM-PCA
    # Add residuals to all LSs (augmented for LiMM-PCA and similar)
    for(i in 1:length(approved)){
      a <- approved[i]
      LS_aug[[effs[a]]] <- LS[[effs[a]]] + residuals
      error[[effs[a]]] <- LS[[effs[a]]] + residuals
    }
    anonam <- rownames(ano)
    # Alphabetically sorted interaction names
    for(i in 1:length(anonam)){
      if(grepl(":", anonam[i], fixed=TRUE)){
        anonam[i] <- paste(sort(strsplit(anonam[i],":")[[1]]), collapse=":")
      }
    }
    dfNum   <- ano[["Df"]]
    dfDenom <- c(rep(ano["Residuals","Df"], length(dfNum)-1),0)
    names(dfNum) <- names(dfDenom) <- anonam # May need to limit to approved?
  } else {
    # Augment errors according to Mixed Model ANOVA and/or LiMM-PCA
    if(!lme4 && !is.logical(REML)){
      ets <- ano$err.terms
      anonam <- rownames(ano$anova)
      # Alphabetically sorted interaction names
      for(i in 1:length(anonam)){
        if(grepl(":", anonam[i], fixed=TRUE)){
          anonam[i] <- paste(sort(strsplit(anonam[i],":")[[1]]), collapse=":")
        }
      }
      names(ets) <- anonam
      dfDenom <- ano$denom.df
      dfNum   <- ano$anova[["Df"]]
      names(dfNum) <- names(dfDenom) <- anonam # May need to limit to approved?
    } else {
      formulaOld <- formula
      formula <- formula(paste0(gsub("(1 | ", "r(", as.character(formula(mod)), fixed = TRUE)[c(2,1,3)], collapse=" "))
      mfPre2 <- mfPre
      mfPre2$REML <- NULL
      # Add unrestricted to mf for use with mixlm
      if(is.null(mfPre2$unrestricted))
        mfPre2$unrestricted <- FALSE
      mod_no_lme4 <- eval(mfPre2, envir = environment())
      ano_no_lme4 <- mixlm::AnovaMix(mod_no_lme4, SStype = SStype)
      formula <- formulaOld
      ets <- ano_no_lme4$err.terms
      no_eff <- rownames(ano_no_lme4$anova)
      # Alphabetically sorted interaction names
      for(i in 1:length(no_eff)){
        if(grepl(":", no_eff[i], fixed=TRUE)){
          no_eff[i] <- paste(sort(strsplit(no_eff[i],":")[[1]]), collapse=":")
        }
      }
      names(ets) <- no_eff
      dfDenom <- ano_no_lme4$denom.df
      dfNum <- ano_no_lme4$anova[["Df"]]
      names(dfNum) <- no_eff
      names(dfDenom) <- no_eff

      # Effective dimensions calculations
      if(use_ED){
        varcor_random <- lapply(models, VarCorr)
        varcor_resid <- sapply(models, sigma)^2
        randEffList <- getME(models[[1]], "flist")
        ED <- list()
        for(i in 1:ncol(Y)){
          wrong.ED <- FALSE
          y <- Y[,i]
          q <- ncol(Z)
          G <- diag(q) * 0
          r0 <- 0
          for(j in 1:length(randEffList)){
            # TODO: nlevR only holds for full coding of random effects
            nlevR <- nlevels(randEffList[[j]])
            G[r0+(1:nlevR), r0+(1:nlevR)] <- diag(nlevR)*
              max(varcor_random[[i]][[names(randEffList[j])]][1], 10^-12)
            r0 <- r0 + nlevR
          }
          R <- diag(N) * varcor_resid[i]
          Rinv <- diag(N) / varcor_resid[i]

          XRX <- t(X) %*% Rinv %*% X
          ZRX <- t(Z) %*% Rinv %*% X
          ZRZ <- t(Z) %*% Rinv %*% Z
          # Catch singular matrix
          try(Ginv <- solve(G), silent = TRUE)
          if(!exists("Ginv")){
            Ginv <- G*Inf
            wrong.ED <- TRUE
            warning("Random effects covariance matrix is singular, using Inf instead")
          }
          ZRZG <- ZRZ + Ginv
          XRY <- t(X) %*% Rinv %*% y
          ZRY <- t(Z) %*% Rinv %*% y

          Q <- solve(rbind(cbind(XRX, t(ZRX)), cbind(ZRX, ZRZG)))
          vecbc <- Q %*% rbind(XRY, ZRY)

          H <- M %*% Q %*% t(M) %*% Rinv
          K <- Q %*% cbind(rbind(XRX, ZRX), rbind(t(ZRX), ZRZ))
          K_diag <- diag(K)[-(1:ncol(X))]
          ED_rand <- numeric(length(randEffList))
          r0 <- 0
          for(j in 1:length(randEffList)){
            nlevR <- nlevels(randEffList[[j]])
            ED_rand[j] <- sum(K_diag[r0+(1:nlevR)])
            r0 <- r0 + nlevR
          }
          if(wrong.ED)
            ED_rand <- rep(NaN, length(ED_rand))
          ED[[i]] <- ED_rand
        }
        EDm <- matrix(unlist(ED), nrow = length(randEffList))
        EDm <- pmax(EDm, 1)
        sortedRandEffList <- names(randEffList)
        # Alphabetically sorted interaction names
        for(i in 1:length(sortedRandEffList)){
          if(grepl(":", sortedRandEffList[i], fixed=TRUE)){
            sortedRandEffList[i] <- paste(sort(strsplit(sortedRandEffList[i],":")[[1]]), collapse=":")
          }
        }
        rownames(EDm) <- sortedRandEffList
        EDall <- matrix(dfNum, nrow = length(dfNum), ncol = ncol(EDm))
        rownames(EDall) <- names(dfNum)
        colnames(EDall) <- colnames(Y)
        EDall["Residuals", ] <- colSums(EDall) + 1 - colSums(EDm)
        EDall[rownames(EDm), ] <- EDm
      }
    }
    for(i in 1:length(approved)){
      a <- approved[i]
      C <- 1

      if(is.numeric(aug_error)){
        alpha <- aug_error
        if(use_ED && any(is.nan(EDall[effs[a],]))){
          warning("Indeterminate effective dimensions, using degrees of freedom instead")
          use_ED <- FALSE
        }
        if(use_ED){
          C <- sqrt(EDall[effs[a],] / EDall[as.numeric(names(ets[[effs[a]]])),]) *
            qf(1-alpha, EDall[effs[a],], EDall[as.numeric(names(ets[[effs[a]]])),])
          C <- matrix(C, nrow=nrow(Y), ncol=ncol(Y), byrow=TRUE)
        } else {
          # LiMM-PCA -> sqrt(dfNum / dfDenom * F(dfNum, dfDenom, 1-alpha))
          #             sqrt(dfNum / dfDenom * pf(1-alpha, dfNum, dfDenom))
          C <- sqrt(dfNum[effs[a]] / dfDenom[effs[a]] * qf(1-alpha, dfNum[effs[a]], dfDenom[effs[a]]))
          C <- matrix(C, nrow(Y), ncol(Y))
        }
      }

      if(as.numeric(names(ets[effs[[a]]][[1]])) == length(effs)+1){
        # Error term is residual
        LS_aug[[effs[a]]] <- LS[[effs[a]]] + residuals * C
        error[[effs[a]]] <- LS[[effs[a]]] + residuals * C
        attr(error[[effs[a]]], 'term') <- "residuals"
      } else {
        # Other error term
        if(exists("no_eff"))
          lsa <- LS[no_eff[as.numeric(names(ets[effs[[a]]][[1]]))]]     # TODO: Check what happens here when adding a numeric effect
        else
          lsa <- LS[as.numeric(names(ets[effs[[a]]][[1]]))]     # TODO: Check what happens here when adding a numeric effect
        # Loop over and sum up (in case of compound errors)
        for(j in 1:length(lsa)){
          error[[effs[a]]] <- error[[effs[a]]] + lsa[[j]] * C
        }
        LS_aug[[effs[a]]] <- LS_aug[[effs[a]]]+error[[effs[a]]]
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

  ########################## Combined effects ##########################
  # Combine effects if the comb() function is used
  eff_combined <- rep(FALSE, length(effs))
  names(eff_combined) <- effs
  remove <- numeric()
  approvedComb <- as.list(approved)
  if(is.list(combined)){
    for(i in 1:length(combined)){
      approved <- c(approved, length(effs)+1)
      approvedAB <- c(approvedAB, length(effs)+1)
      approvedComb[[length(approvedComb)+1]] <- numeric()
      combName <- paste(combined[[i]], collapse="+")
      names(approved)[length(approved)] <- names(approvedAB)[length(approvedAB)] <- combName
      LS[[approved[length(approved)]]] <- 0*LS[[effs[[approved[1]]]]]
      names(LS)[approved[length(approved)]] <- combName
      effs[approved[length(approved)]] <- combName
      eff_combined[approved[length(approved)]] <- TRUE
      names(eff_combined)[approved[length(approved)]] <- combName
      error[[approved[length(approved)]]] <- 0*residuals
      names(error)[approved[length(approved)]] <- combName
      ssq[approved[length(approved)]] <- 0
      names(ssq)[approved[length(approved)]] <- combName
      maxDir[approved[length(approved)]] <- 0
      for(dis in combined[[i]]){
        if(grepl(":", dis, fixed=TRUE)) # Reorder interactions to alphabetical order
          dis <- paste(sort(strsplit(dis,":")[[1]]), collapse=":")
        if(any(!(dis %in% names(approvedAB))))
          stop("Cannot combine a continuous effect with a categorical factor.")
        remove <- c(remove, which(effsAB==dis))
        # Remove from approved the entry named dis
        LS[[approved[length(approved)]]] <- LS[[approved[length(approved)]]] + LS[[effs[approvedAB[names(approvedAB)==dis]]]]
        ssq[approved[length(approved)]] <- ssq[approved[length(approved)]] + ssq[effs[approvedAB[names(approvedAB)==dis]]]
        maxDir[approved[length(approved)]] <- max(maxDir[approved[length(approved)]], maxDir[approvedAB[names(approvedAB)==dis]])
        # Accumulate combined effect numbers
        approvedComb[[length(approvedComb)]] <- c(approvedComb[[length(approvedComb)]], approvedAB[names(approvedAB)==dis])
        approved <- approved[names(approvedAB)!=dis]
        approvedAB <- approvedAB[names(approvedAB)!=dis]
      }
      # Assume residual error for combined effect?
      error[[approved[length(approved)]]] <- LS[[approved[length(approved)]]] + residuals
      #      error[[approved[length(approved)]]] <- error[[approved[length(approved)]]] + error[[effs[approvedAB[names(approvedAB)==dis]]]]
      LS_aug[[approved[length(approved)]]] <- error[[approved[length(approved)]]]
      names(LS_aug)[approved[length(approved)]] <- combName
    }
  } else {
    names(approvedComb) <- names(approvedAB)
  }
  if(names(ssq)[length(ssq)] == "Residuals")
    ssq <- c(ssq[setdiff(1:(length(ssq)-1), remove)],ssq_residual)
  else {
    ssq <- c(ssq[setdiff(1:(length(ssq)), remove)],ssq_residual)
  }
  names(ssq)[length(ssq)] <- "Residuals"
  eff_combined <- eff_combined[approved]

  ########################## Permutation ##########################
  # Permutation testing
  if(!(is.logical(permute) && !permute)){
    # Default to 1000 permutations
    if(is.logical(permute)){
      permute <- 1000
      if(interactive())
        message("Defaulting to 1000 permutations\n")
    }
    if(perm.type[1] == "approximate"){
      ssqa <- pvals <- numeric(length(approved))
      ssqaperm <- lapply(1:length(approved),function(i)"NA")
      names(ssqa) <- names(ssqaperm) <- names(pvals) <- effs[approved]
      for(i in 1:length(approved)){
        a <- approvedAB[i]
        perms <- numeric(permute)
        # Subset of design matrix for effect a
        D <- M[, assign%in%approvedComb[[names(a)]], drop=FALSE]
        DD <- D %*% pracma::pinv(D)
        # Base ssq
        ssqa[effs[a]] <- norm(DD %*% LS_aug[[effs[a]]], "F")^2
        pb <- progress_bar$new(total = permute, format = paste0("  Permuting ", effs[a], " (", i,"/",length(approved),") [:bar] :percent (:eta)"))
        # Permuted ssqs
        for(perm in 1:permute){
          perms[perm] <- norm(DD %*% LS_aug[[effs[a]]][sample(N),], "F")^2
          pb$tick()
        }
        ssqaperm[[effs[a]]] <- perms
        pvals[effs[a]] <- sum(perms > ssqa[effs[a]])/(permute)
      }
    }
    if(perm.type[1] == "exact"){
      # lme4 vs LS
      if(!lme4 && !is.logical(REML)){
        stop("Exact permutations in lme4 models not implemented yet!")
      }
      # Loop over effects
      for(i in 1:length(approved)){
        stop("Exact permutations not implemented yet!")
        a <- approvedAB[i]
        perms <- numeric(permute)
        # Subset of design matrix for effect a
        #D <- M[, assign%in%approvedComb[[names(a)]], drop=FALSE]
        #DD <- D %*% pracma::pinv(D)
      # Find permissible units
        # Check balance, warning
      # Permute
      }
    }
  }

  ########################## SCA ##########################
  # SCAs
  scores <- loadings <- projected <- singulars <- list()
  for(i in approved){
    maxDiri <- min(Rank(LS[[effs[i]]]),maxDir[i])
    if(pca.in != 0)
      maxDiri <- min(maxDiri, pca.in)
    if(add_error)
      maxDiri <- min(N-1, p)
    if(maxDiri == 0)
      stop(paste0("Effect '", effs[i], "' has no estimable levels"))
    pcai <- .pca(LS[[effs[i]]], ncomp=maxDiri, proj=error[[effs[i]]])
    scores[[effs[i]]] <- pcai$scores
    loadings[[effs[i]]] <- pcai$loadings
    projected[[effs[i]]] <- pcai$projected
    singulars[[effs[i]]] <- pcai$singulars

    if(pca.in!=0){ # Transform back if PCA on Y has been performed
      loadings[[effs[i]]] <- pca$loadings[,1:pca.in,drop=FALSE] %*% loadings[[effs[i]]]
      dimnames(loadings[[effs[i]]]) <- list(colnames(Yorig), paste("Comp", 1:maxDiri, sep=" "))
    }
  }
  # SCA of residuals
  maxDirRes <- min(N-1,p)
  if(pca.in != 0)
    maxDirRes <- min(maxDirRes, pca.in)
  pcaRes <- .pca(residuals, ncomp=maxDirRes)
  scores[["Residuals"]] <- pcaRes$scores
  loadings[["Residuals"]] <- pcaRes$loadings
  projected[["Residuals"]] <- pcaRes$projected
  singulars[["Residuals"]] <- pcaRes$singulars

  # Create model.frame object
  model <- model.frame(mod)
  model[[1]] <- Yorig

  ########################## Return ##########################
  obj <- list(scores=scores, loadings=loadings, projected=projected, singulars=singulars,
              LS=LS, effects=effects, coefficients=coefs, Y=Yorig, X=M, residuals=residuals,
              error=error, eff_combined=eff_combined, SStype=SStype, contrasts=contrasts, unrestricted=unrestricted,
              ssq=ssq, ssqY=ssqY, explvar=ssq/ssqY, models=models, anovas=anovas, model.frame=modFra,
              call=match.call(), fit.type=fit.type, add_error=add_error, dfNum=dfNum, dfDenom=dfDenom,
              model = model)
  if(pca.in!=0){
    obj$Ypca <- list(pca=pca, ncomp=pca.in)
  }
  if(permute || is.numeric(permute)){
    obj$permute <- list(ssqa=ssqa, ssqaperm=ssqaperm, pvalues=pvals, permutations=permute)
  }
  if(use_ED)
    obj$ED <- EDall
  if(exists("ets")){
    denoms <- unlist(lapply(ets, function(e){ifelse(is.na(e),NA,as.numeric(names(e)))}))
    names(denoms) <- names(ets)
    obj$denoms <- denoms
  } else {
    denoms <- c(rep(which(names(dfDenom) == "Residuals"), length(dfDenom)-1),NA)
    names(denoms) <- names(dfDenom)
    obj$denoms <- denoms
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

# Not used
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
cparse <- function (f, REML = NULL) {
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
  if(is.logical(REML)){
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
