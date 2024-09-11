#' @title ASCA Plot Methods
#' @name asca_plots
#' @aliases asca_plots scoreplot.asca loadingplot.asca
#'
#' @description Various plotting procedures for \code{\link{asca}} objects.
#'
#' @details Usage of the functions are shown using generics in the examples in \code{\link{asca}}.
#' Plot routines are available as
#' \code{scoreplot.asca} and \code{loadingplot.asca}.
#'
#' @param object \code{asca} object.
#' @param factor \code{integer/character} for selecting a model factor. If factor <= 0 or "global",
#' the PCA of the input is used (negativ factor to include factor level colouring with global PCA).
#' @param comps \code{integer} vector of selected components.
#' @param pch.scores \code{integer} plotting symbol.
#' @param pch.projections \code{integer} plotting symbol.
#' @param gr.col \code{integer} vector of colours for groups.
#' @param ellipsoids \code{character} "confidence" or "data" ellipsoids for balanced fixed effect models.
#' @param confidence \code{numeric} vector of ellipsoid confidences, default = c(0.4, 0.68, 0.95).
#' @param xlim \code{numeric} x limits.
#' @param ylim \code{numeric} y limits.
#' @param xlab \code{character} x label.
#' @param ylab \code{character} y label.
#' @param legendpos \code{character} position of legend.
#' @param ... additional arguments to underlying methods.
#'
#' @return The plotting routines have no return.
#'
#' @references
#' * Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J., and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA): A new tool for analyzing designed metabolomics data. Bioinformatics, 21(13), 3043–3048.
#' * Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence ellipsoids for ASCA models based on multivariate regression theory. Journal of Chemometrics, 32(e2990), 1–13.
#' * Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and linear mixed models to analyse high-dimensional designed data. Journal of Chemometrics, 34(6), e3232.
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for computation and extraction of results are found in \code{\link{asca_results}}.
#'
#' @export
loadingplot.asca <- function(object, factor = 1, comps = 1:2, ...){
  if(inherits(object,"pcanova")){
    # Remove too high component numbers
    comps <- comps[comps <= length(object$anovas)]
    if(length(comps) == 0)
      stop("No components to plot")
    if(factor > 0)
      factor <- -factor
  }
  if((inherits(object, "scores") && ncol(object) == 1 ) ||
     (inherits(object, "multiblock") && length(object$Xvar) == 1)){ # Check for single component in model
    comps <- comps[1]
    nComps <- length(comps)
  }
  # Check if input PCA should be used instead of PCAs of factor LS matrices.
  global <- FALSE
  if(factor < 1 || factor == "global"){
    if(!is.null(object$Ypca)){
      loads <- object$Ypca$pca$loadings
      global <- TRUE
      factor <- abs(factor)
    } else {
      warning("No global PCA available, using first factor")
      factor <- 1
      global <- FALSE
    }
  }
  if(!global){
    loads <- loadings(object=object, factor=factor)
  }
  plot(loads, comps=comps, ...)
}

#' @rdname asca_plots
#' @export
scoreplot.asca <- function(object, factor = 1, comps = 1:2, within_level = "all",
                           pch.scores = 19, pch.projections = 1,
                           gr.col = NULL, projections = TRUE,
                           spider = FALSE, ellipsoids, confidence,
                           xlim,ylim, xlab,ylab, legendpos, ...){

  if(factor == "within" || factor == "Within")
    factor <- "Residuals"
  if(is.null(gr.col)){
    if(inherits(object,"msca"))
      gr.col <- 1:nlevels(object$effects[[1]])
    else
      gr.col <- 1:nlevels(object$effects[[factor]])
  }
  if(inherits(object,"pcanova")){
    # Remove too high component numbers
    comps <- comps[comps <= length(object$anovas)]
    if(length(comps) == 0)
      stop("No components to plot")
    if(factor > 0)
      factor <- -factor
    object$add_error <- TRUE
    if(factor == 0)
      object$add_error <- FALSE
  }
  # Check if input PCA should be used instead of PCAs of factor LS matrices.
  global <- FALSE
  if(factor < 1 || factor == "global"){
    if(!is.null(object$Ypca)){
      scors <- projs <- object$Ypca$pca$scores
      global <- TRUE
      factor <- abs(factor)
    } else {
      warning("No global PCA available, using first factor")
      factor <- 1
      global <- FALSE
    }
  }

  nobj <- nrow(object$Y)
  nlev <- 0
  if(!(factor==0) && !global && !factor == "Residuals" && !factor == length(object$scores)){
    # Number of levels in current factor
    nlev  <- nlevels(object$effects[[factor]])
    # Remove redundant levels
    comps <- comps[comps <= nlev-1]
    # Set gr.col if missing
    if(missing(gr.col)){
      gr.col <- adjustcolor(rep(palette(), max(1, nlev%/%8+1))[1:nlev], alpha.f = 0.7)
    }
  }
  if(!global){
    # Get scores and projections
    projs <- scors <- scores(object=object, factor=factor)
    if(!inherits(object,"pcanova")){
      if(factor != "Residuals" && factor != length(object$scores))
        projs <- projections(object=object, factor=factor) #+ scors
    }
  }

  # Extract scores for within factor if MSCA with single level within factor
  if(inherits(object, "msca") &&
     (factor == length(object$scores) || factor == "Residuals") &&
     within_level != "all"){
    projs <- scors <- object$scores.within[[within_level]]
  }

  # Find limits of plotting
  if(missing(xlim))
    xlim <- c(min(min(scors[,comps[1]]), min(projs[,comps[1]])),
              max(max(scors[,comps[1]]), max(projs[,comps[1]])))
  if(missing(ylim))
    if(length(comps)>1)
      ylim <- c(min(min(scors[,comps[2]]), min(projs[,comps[2]])),
                max(max(scors[,comps[2]]), max(projs[,comps[2]])))
  else
    ylim <- c(0.5, nlev+0.5)
  # Generate labels
  evar <- attr(scors, 'explvar')
  if(missing(xlab))
    xlab <- paste0("Comp ", comps[1], " (",format(evar[comps[1]], digits = 2, trim = TRUE), " %)")
  if(missing(ylab))
    if(length(comps)>1)
      ylab <- paste0("Comp ", comps[2], " (",format(evar[comps[2]], digits = 2, trim = TRUE), " %)")
  else
    ylab <- 'Level'


  # Scatter plot
  if(length(comps)>1){
    if(object$add_error || inherits(object, "msca")) # Skip plotting of scores if error is added (APCA)
      scoreplot(scors, comps=comps, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch.projections, col="white", ...)
    else
      scoreplot(scors, comps=comps, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, pch=pch.scores, ...)
    # Add projections
    if(factor!=0 && !global && !factor == "Residuals" && !factor == length(object$scores))
      for(i in 1:nlev){
        lev <- levels(object$effects[[factor]])[i]
        if(!object$add_error) # Skip plotting of scores if error is added (APCA)
          points(scors[object$effects[[factor]] == lev, comps], pch=pch.scores, col=gr.col[i])
        # Backprojections
        if(projections && !(factor==0)) # Skip projections if global PCA is used
          points(projs[object$effects[[factor]] == lev, comps], pch=pch.projections, col=gr.col[i])
        # Spider plot
        if(spider && !(factor==0)){
          score_j <- scors[object$effects[[factor]] == lev, comps]
          projs_j <- projs[object$effects[[factor]] == lev, comps]
          for(j in 1:nrow(projs_j)){
            lines(c(score_j[j,1], projs_j[j,1]), c(score_j[j,2], projs_j[j,2]), col=gr.col[i])
          }
        }
      }
    # MSCA scoreplot for within factor
    if(inherits(object, "msca") && (factor == length(object$scores) || factor == "Residuals")){
      scors1 <- scores(object=object, factor=1)
      nlev1 <- nlevels(object$effects[[1]])
      if(is.numeric(within_level))
        w_levs <- within_level
      if(within_level == "all")
        w_levs <- 1:nlev1
      for(i in w_levs){
        lev1 <- levels(object$effects[[1]])[i]
        if(projections && !(factor==0) && within_level != "all") # Skip projections if global PCA is used
          points(projs[, comps], pch=pch.projections, col=gr.col[1])
        if(projections && !(factor==0) && within_level == "all") # Skip projections if global PCA is used
          points(projs[object$effects[[1]] == lev1, comps], pch=pch.projections, col=gr.col[i])
      }
    }
    if(!missing(legendpos))
      if(factor!=0)
        legend(legendpos, legend = levels(object$effects[[factor]]), col=gr.col, pch=pch.scores)

    if(!missing(ellipsoids)){
      if(object$eff_combined[factor])
        stop("Ellipsoids not defined for combined effects")
      if(missing(confidence))
        confidence <- c(0.4,0.68,0.95)
      if(ellipsoids == "data"){
        dataEllipse(projs[,comps], groups = object$effects[[factor]], levels=confidence, add=TRUE, plot.points=FALSE, col=rep("black", length(gr.col)), lwd=1, group.labels="", center.pch=FALSE, lty=1)
        dataEllipse(projs[,comps], groups = object$effects[[factor]], levels=confidence, add=TRUE, plot.points=FALSE, col=gr.col, lwd=1, group.labels="", center.pch=FALSE, lty=2)
      }
      if(ellipsoids == "confidence" || ellipsoids == "conf" || ellipsoids == "model"){
        if(inherits(object, "apca"))
          stop("Confidence ellipsoids not implemented for APCA")
        # Covariance matrix
        sigma <- crossprod(object$error[[factor]]-object$LS[[factor]])/nobj
        #        sigma <- crossprod(object$residuals)/nobj
        L <- object$loadings[[factor]][,comps]
        # Transformed covariance matrix
        LSL <- crossprod(L,sigma) %*% L * nlev / (nobj-nlev)
        if(!isSymmetric(LSL))
          LSL <- (LSL + t(LSL))/2 # Force symmetry
        # Scaling by confidence
        cx <- list()
        for(c in 1:length(confidence))
          cx[[c]] <- sqrt((nobj-nlev)*2 / (nobj-nlev-2+1) * qf(confidence[c], 2, nobj-nlev-2+1))
        for(i in 1:nlev){
          lev <- levels(object$effects[[factor]])[i]
          for(c in 1:length(confidence)){
            ellipse(colMeans(scors[object$effects[[factor]]==lev,comps]), LSL, cx[[c]], lwd=1, col="black", center.pch=FALSE)
            ellipse(colMeans(scors[object$effects[[factor]]==lev,comps]), LSL, cx[[c]], lwd=1, col=gr.col[i], lty=2, center.pch=FALSE)
          }
        }
      }
    }
  } else { # Line plot
    plot(scors[,comps], as.numeric(object$effects[[factor]]), xlim=xlim,
         ylim=ylim, xlab=xlab, ylab=ylab, axes = FALSE)
    axis(1)
    axis(2, at=1:nlev, labels = levels(object$effects[[factor]]))
    box()
    for(i in 1:nlev){
      lev <- levels(object$effects[[factor]])[i]
      if(!object$add_error) # Skip plotting of scores if error is added (APCA/LiMM-PCA)
        points(scors[object$effects[[factor]] == lev, comps], rep(i,sum(as.numeric(object$effects[[factor]]) == i)), pch=pch.scores, col=gr.col[i])
      if(!(factor==0)) # Skip projections if global PCA is used
        points(projs[object$effects[[factor]] == lev, comps], rep(i,sum(as.numeric(object$effects[[factor]]) == i)), pch=pch.projections, col=gr.col[i])
    }
    if(!missing(ellipsoids)){
      if(missing(confidence))
        confidence <- c(0.4,0.68,0.95)
      if(ellipsoids == "confidence" || ellipsoids == "conf" || ellipsoids == "model"){
        sigma <- crossprod(object$error[[factor]]-object$LS[[factor]])/nobj
        # sigma <- crossprod(object$residuals)/nobj
        L <- object$loadings[[factor]][,comps]
        # Transformed covariance matrix
        LSL <- sqrt(crossprod(L,sigma) %*% L * nlev / (nobj-nlev))
        # Scaling by confidence
        cx <- list()
        for(c in 1:length(confidence))
          cx[[c]] <- sqrt((nobj-nlev)*1 / (nobj-nlev-1+1) * qf(confidence[c], 1, nobj-nlev-1+1))
        for(i in 1:nlev){
          lev <- levels(object$effects[[factor]])[i]
          for(c in 1:length(confidence)){
            lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)+c(LSL)*cx[[c]], i+c(-0.2,0.2), col=gr.col[i])
            lines(mean(scors[object$effects[[factor]] == lev,comps])*c(1,1)-c(LSL)*cx[[c]], i+c(-0.2,0.2), col=gr.col[i])
          }
        }
      }
    }
  }
}

#' @rdname asca_plots
#' @export
permutationplot <- function(object, factor = 1, xlim, xlab = "SSQ", main, ...){
  if(is.null(object$permute))
    stop("Re-run model with permutation testing enabled")
  if(missing(xlim))
    xlim <- range(c(object$permute$ssqaperm[[factor]], object$ssq[factor]))
  if(missing(main)){
    if(is.numeric(factor))
      main <- paste0("Permutation of ", names(object$ssq)[factor], " effect")
    else
      main <- paste0("Permutation of ", factor, " effect")
  }
  hist(object$permute$ssqaperm[[factor]], xlim=xlim, xlab=xlab, main=main, ...)
  abline(v = object$ssq[factor], col=2, lwd=2, ...)
}

#' @title Timeplot for Combined Effects
#' @name timeplot
#'
#' @param object \code{asca} object.
#' @param factor \code{integer/character} main factor.
#' @param time \code{integer/character} time factor.
#' @param comb \code{integer/character} combined effect factor.
#' @param comp \code{integer} component number.
#' @param ylim \code{numeric} y limits.
#' @param x_time \code{logical} use time levels as non-equispaced x axis (default = FALSE).
#' @param xlab \code{character} x label.
#' @param ylab \code{character} y label.
#' @param lwd \code{numeric} line width.
#' @param ... additional arguments to \code{plot}.
#'
#' @return Nothing
#' @export
#'
#' @examples
#' data("caldana")
#' mod.comb <- asca(compounds ~ time + comb(light + light:time), data=caldana)
#'
#' # Default time axis
#' timeplot(mod.comb, factor="light", time="time", comb=2)
#'
#' # Non-equispaced time axis (using time levels)
#' timeplot(mod.comb, factor="light", time="time", comb=2, x_time=TRUE)
#'
#' # Second component
#' timeplot(mod.comb, factor="light", time="time", comb=2, comp=2, x_time=TRUE)
timeplot <- function(object, factor, time, comb, comp=1, ylim, x_time = FALSE,
                     xlab = time, ylab = paste0("Score ",comp), lwd = 2, ...){
  # Extract levels for time and factor
  time_fac <- as.numeric(as.character(object$model.frame[[time]]))
  if(x_time)
    x <- as.numeric(object$model.frame[[time]])
  else
    x <- time_fac
  factor_fac <- object$model.frame[[factor]]
  comb_scores <- object$scores[[comb]][ ,comp]
  if(missing(ylim))
    ylim <- c(min(comb_scores),max(comb_scores))
  plot(x[factor_fac==levels(factor_fac)[1]],
       comb_scores[factor_fac==levels(factor_fac)[1]], type="l",
       ylab = ylab, xlab = xlab,
       ylim = ylim, lwd=lwd,
       panel.first = {abline(h=0, lty=2, col="gray")}, axes=!x_time, ...)
  for(i in 2:length(levels(factor_fac))){
    lines(x[factor_fac==levels(factor_fac)[i]],
          comb_scores[factor_fac==levels(factor_fac)[i]],
          type="l", col=i, lty=i, lwd=lwd)
  }
  if(x_time){
    axis(1, at=1:length(levels(object$model.frame[[time]])), labels = levels(object$model.frame[[time]]))
    axis(2)
    box()
  }
}
