#' @title ASCA Result Methods
#' @name asca_results
#' @aliases asca_results print.asca summary.asca projections projections.asca print.summary.asca loadings.asca scores.asca
#'
#' @description Standard result computation and extraction functions for ASCA (\code{\link{asca}}).
#'
#' @details Usage of the functions are shown using generics in the examples in \code{\link{asca}}.
#' Explained variances are available (block-wise and global) through \code{blockexpl} and \code{print.rosaexpl}.
#' Object printing and summary are available through:
#' \code{print.asca} and \code{summary.asca}.
#' Scores and loadings have their own extensions of \code{scores()} and \code{loadings()} through
#' \code{scores.asca} and \code{loadings.asca}. Special to ASCA is that scores are on a
#' factor level basis, while back-projected samples have their own function in \code{projections.asca}.
#'
#' @param object \code{asca} object.
#' @param x \code{asca} object.
#' @param factor \code{integer/character} for selecting a model factor.
#' @param digits \code{integer} number of digits for printing.
#' @param ... additional arguments to underlying methods.
#'
#' @return Returns depend on method used, e.g. \code{projections.asca} returns projected samples,
#' \code{scores.asca} return scores, while print and summary methods return the object invisibly.
#'
#' @references
#' * Smilde, A., Jansen, J., Hoefsloot, H., Lamers,R., Van Der Greef, J., and Timmerman, M.(2005). ANOVA-Simultaneous Component Analysis (ASCA): A new tool for analyzing designed metabolomics data. Bioinformatics, 21(13), 3043–3048.
#' * Liland, K.H., Smilde, A., Marini, F., and Næs,T. (2018). Confidence ellipsoids for ASCA models based on multivariate regression theory. Journal of Chemometrics, 32(e2990), 1–13.
#' * Martin, M. and Govaerts, B. (2020). LiMM-PCA: Combining ASCA+ and linear mixed models to analyse high-dimensional designed data. Journal of Chemometrics, 34(6), e3232.
#'
#' @seealso Overviews of available methods, \code{\link{multiblock}}, and methods organised by main structure: \code{\link{basic}}, \code{\link{unsupervised}}, \code{\link{asca}}, \code{\link{supervised}} and \code{\link{complex}}.
#' Common functions for plotting are found in \code{\link{asca_plots}}.
#'
#' @export
print.asca <- function(x, ...){
  mod <- "Anova Simultaneous Component Analysis"
  if(inherits(x, "apca"))
    mod <- "Anova Principal Component Analysis"
  if(inherits(x, "limmpca"))
    mod <- "LiMM-PCA"
  if(inherits(x, "msca"))
    mod <- "Multilevel Simultaneous Component Analysis"
  cat(paste0(mod, " fitted using"), x$fit.type)
  cat("\nCall:\n", deparse(x$call), "\n", sep = "")
  invisible(x)
}

#' @rdname asca_results
#' @export
summary.asca <- function(object, extended=TRUE, df=FALSE, ...){
  dat <- data.frame(SSQ=object$ssq, "Expl.var"=object$explvar*100)
  colnames(dat) <- c("Sum.Sq.", "Expl.var.(%)")
#  dat <- dat[-nrow(dat),,drop=FALSE]
  if(!is.null(object$permute)){
    pvals <- object$permute$pvalues
    pvals[pvals==0] <- 1/object$permute$permutations
    pv <- rep(NA,nrow(dat))
    names(pv) <- rownames(dat)
    pv[names(pvals)] <- pvals
    dat <- cbind(dat, "P-value"=pv)
  }
  mod <- "Anova Simultaneous Component Analysis"
  if(inherits(object, "apca"))
    mod <- "Anova Principal Component Analysis"
  if(inherits(object, "limmpca"))
    mod <- "LiMM-PCA"
  if(inherits(object, "msca")){
    mod <- "Multilevel Simultaneous Component Analysis"
    rownames(dat) <- c("Between", "Within")
  }
  x <- list(dat=dat, model=mod, fit.type=object$fit.type)
  if(extended){
    LS_REML <- "least squares"
    if(!inherits(object$models[[1]],"lm"))
      LS_REML <- ifelse(getME(object$models[[1]],"is_REML"), "REML", "ML")
    ss <- c("I","II","III")
    x$info <- paste0("SS type ", ss[object$SStype], ", ", object$coding, " coding, ",
                     ifelse(object$unrestricted, "unrestricted","restricted"), " model",
                     ", ", LS_REML, " estimation")
    if(!is.null(object$permute))
      x$info <- paste0(x$info, ", ", object$permute$permutations, " permutations")
  }
  if(df){
    x$dat <- cbind(x$dat, "df"=object$dfNum, "df.denom"=object$dfDenom, "err.term"=object$denom)
  }
  class(x) <- c('summary.asca')
  x
}

#' @rdname asca_results
#' @export
print.summary.asca <- function(x, digits=2, ...){
  cat(x$mod, "fitted using", x$fit.type, "\n")
  if(!is.null(x$info))
    cat("-", x$info, "\n")
  print(round(x$dat, digits))
  invisible(x$dat)
}

#' @rdname asca_results
#' @export
loadings.asca <- function(object, factor = 1, ...){
  loads <- object$loadings[[factor]]
  class(loads) <- "loadings"
  return(loads)
}

#' @rdname asca_results
#' @export
scores.asca <- function(object, factor = 1, ...){
  scors <- object$scores[[factor]]
  class(scors) <- "scores"
  return(scors)
}

#' @rdname asca_results
#' @export
projections <- function (object, ...) {
  UseMethod("projections", object)
}

#' @rdname asca_results
#' @export
projections.asca <- function(object, factor = 1, ...){
  projs <- object$projected[[factor]]
  class(projs) <- "projs"
  return(projs)
}
