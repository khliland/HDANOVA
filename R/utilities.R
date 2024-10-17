# PCA for internal use
.pca <- function(X, ncomp, scale=FALSE, proj=NULL, ...){
  X <- as.matrix(unclass(X))
  if(!inherits(X,'matrix'))
    stop("'X' must be a matrix")
  if(missing(ncomp))
    ncomp <- min(c(nrow(X)-1, ncol(X)))
  y   <- structure(matrix(rnorm(nrow(X)), ncol=1), rownames = rownames(X))
  dat <- data.frame(y=y, X = I(X))
  mod <- pls::pcr(y ~ X, ncomp = ncomp, data = dat, scale = scale, ...)
  mod$explvar <- mod$Xvar/mod$Xtotvar*100
#  PCR <- pls::pcr(y ~ X, ncomp = ncomp, data = dat, scale = scale, ...)
#  mod <- list(loadings=PCR$loadings, scores=PCR$scores, Xmeans=PCR$Xmeans, explvar=PCR$Xvar/PCR$Xtotvar*100, PCA = PCR)
  attr(mod$loadings, 'explvar') <- attr(mod$scores, 'explvar') <- mod$explvar
  if(!is.null(proj)){
    mod$projected <- (proj-rep(mod$Xmeans, each=nrow(proj))) %*% mod$loadings
  }
  mod$singulars <- sqrt(mod$Xvar)
  mod$call <- match.call()
  class(mod) <- c('pca')
  return(mod)
}

# Cross-validated PCA (Eigenvector approach)
.PCAcv <- function(X, ncomp){

  X <- as.matrix(X)
  class(X) <- "matrix"
  N <- dim(X)

  # Center X
  X <- X - rep(colMeans(X), each = N[1])

  # Set ncomp
  max.comp <- min(N[1]-1,N[2])
  if(missing(ncomp)){
    ncomp <- max.comp
  } else {
    ncomp <- min(ncomp, max.comp)
  }

  # Prepare storage
  Xhat <- array(0, dim = c(N[1],N[2],ncomp))

  # Select svd version
  if(ncomp < max.comp){
    sv <- svds
  } else {
    sv <- function(X, k, nv, nu)
      svd(X, nv = nv, nu = nu)
  }

  # Cross-validation (leave-one-out)
  pb <- progress_bar$new(total = N[1], format = "  [:bar] :percent (:eta)")
  for(i in 1:N[1]){
    Xi  <- X[-i, , drop = FALSE]
    Pi  <- sv(Xi, k = ncomp,  nv = ncomp, nu = 0)$v
    Xii <- matrix(rep(X[i,], N[2]), N[2], N[2], byrow = TRUE)
    diag(Xii) <- 0

    # Magic to avoid information bleed
    PiP  <- apply(Pi^2, 1, cumsum)
    PiP1 <- t(PiP/(1-PiP)+1)
    PihP <- t(Pi*(Xii%*%Pi))
    for(j in 1:N[2]){
      PP <- PihP[,j, drop = FALSE] %*% PiP1[j,, drop = FALSE]
      PP[lower.tri(PP)] <- 0
      Xhat[i,j, ] <- colSums(PP)
    }
    pb$tick()
  }

  error <- numeric(ncomp)
  for(i in 1:ncomp){
    error[i] <- sum((X-Xhat[,,i])^2)
  }
  error
}

#' Extracting the Extended Model Frame from a Formula or Fit
#'
#' @description
#' This function attempts to apply \code{\link{model.frame}} and extend the
#' result with columns of interactions.
#'
#' @param formula a model formula or terms object or an R object.
#' @param data a data.frame, list or environment (see \code{\link{model.frame}}).
#' @param ... further arguments to pass to \code{\link{model.frame}}.
#' @param sep separator in contraction of names for interactions (default = ".").
#'
#' @return A \code{\link{data.frame}} that includes everything a \code{\link{model.frame}}
#' does plus interaction terms.
#' @export
#'
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{asca_fit}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#' @examples
#' dat <- data.frame(Y = c(1,2,3,4,5,6),
#'                   X = factor(LETTERS[c(1,1,2,2,3,3)]),
#'                   W = factor(letters[c(1,2,1,2,1,2)]))
#' extended.model.frame(Y ~ X*W, dat)
extended.model.frame <- function(formula, data, ..., sep = "."){
  # Model frame and terms factors (remove response row if presents)
  fac <- attr(terms(mf <- model.frame(formula, data, ...)), "factors")
  fac <- fac[rowSums(fac)>0,, drop=FALSE]
  # Create all included interactions
  int <- lapply(1:ncol(fac), function(i) interaction(data[rownames(fac)[fac[,i]==1]], sep = sep))
  # Convert to data.frame and add interactions to original model.frame
  names(int) <- colnames(fac)
  modFra <- as.data.frame(int, check.names = FALSE)
  for(i in setdiff(colnames(modFra), colnames(mf)))
    mf[[i]] <- modFra[[i]]
  mf
}


#' Update a Model without Factor
#' @description
#' Perform a model update while removing a chosen factor. Hierarchical
#' corresponds to type "II" sum-of-squares, i.e., obeying marginality,
#' while non-hierarchical corresponds to type "III" sum-of-squares.
#'
#' @return An updated model object is returned. If the supplied model is of
#' type \code{lmerMod} and no random effects are left, the model is
#' automatically converted to a linear model before updating.
#'
#' @param model \code{model} object to update.
#' @param fac \code{character} factor to remove.
#' @param hierarchical \code{logical} obey hierarchy when removing factor (default = TRUE).
#' @importFrom lme4 findbars
#' @export
update_without_factor <- function(model, fac, hierarchical = TRUE){
  # Extract formula from model
  form <- formula(model)
  if(hierarchical){
    # Remove factor and all interactions containing factor from formula
    facs <- attr(terms(form),"term.labels")
    for(i in grep(fac, facs))
      form <- update(form, paste0(". ~ . - (", facs[i],")"))
  } else {
    # Remove factor from formula
    form <- update(form, paste0(". ~ . - ", fac, " - (1|", fac, ")"))
  }

  if(inherits(model,"lmerMod") && is.null(findbars(form))){
    # If model is of lme4 type and no random effects are left, remodel as lm() instead
    ca <- getCall(model)
    ca[[1]] <- as.name("lm")
    env <- environment(ca[[2]])
    ca[[2]] <- form
    environment(ca[[2]]) <- env
    eval(ca)
  } else {
    # Update model with new formula
    update(model, form)
  }
}


#' @title Block-wise indexable data.frame
#'
#' @description This is a convenience function for making \code{data.frame}s that are easily
#' indexed on a block-wise basis.
#'
#' @param X Either a single \code{data.frame} to index or a \code{list} of matrices/data.frames
#' @param block_inds Named \code{list} of indexes if \code{X} is a single \code{data.frame}, otherwise \code{NULL}.
#' @param to.matrix \code{logical} indicating if input list elements should be converted to matrices.
#'
#' @return A \code{data.frame} which can be indexed block-wise.
#' @seealso Main methods: \code{\link{asca}}, \code{\link{apca}}, \code{\link{limmpca}}, \code{\link{msca}}, \code{\link{pcanova}}, \code{\link{prc}} and \code{\link{permanova}}.
#' Workhorse function underpinning most methods: \code{\link{asca_fit}}.
#' Extraction of results and plotting: \code{\link{asca_results}}, \code{\link{asca_plots}}, \code{\link{pcanova_results}} and \code{\link{pcanova_plots}}
#' @examples
#' # Random data
#' M <- matrix(rnorm(200), nrow = 10)
#' # .. with dimnames
#' dimnames(M) <- list(LETTERS[1:10], as.character(1:20))
#'
#' # A named list for indexing
#' inds <- list(B1 = 1:10, B2 = 11:20)
#'
#' X <- block.data.frame(M, inds)
#' str(X)
#'
#' @export
block.data.frame <- function(X, block_inds = NULL, to.matrix = TRUE){
  if(!is.null(block_inds)){
    # Use indices to convert matrix/data.frame into list of matrices/data.frames
    if(is.null(names(block_inds)))
      warning("When 'block_inds' is supplied, it should be a named list with indices/names of variables associated with blocks.")
    Z <- lapply(block_inds, function(i)X[,i,drop=FALSE])
    X <- lapply(Z, function(z){rownames(z) <- rownames(X);z})
  }
  # Enclose blocks in "as.is"
  if(to.matrix)
    X <- lapply(X, function(x)I(as.matrix(x)))
  else
    X <- lapply(X, function(x)I(x))
  # Return as data.frame
  X <- do.call(data.frame, X)
  X <- lapply(X, function(x){rownames(x) <- rownames(X);x})
  X <- data.frame(X)
  return(X)
}


#' Dummy-coding of a single vector
#'
#' @description Flexible dummy-coding allowing for all R's built-in types of contrasts
#' and optional dropping of a factor level to reduce rank defficiency probability.
#'
#' @param Y \code{vector} to dummy code.
#' @param contrast Contrast type, default = "contr.sum".
#' @param drop \code{logical} indicating if one level should be dropped (default = TRUE).
#'
#' @return \code{matrix} made by dummy-coding the input vector.
#'
#' @examples
#' vec <- c("a","a","b","b","c","c")
#' dummycode(vec)
#' @export
dummycode <- function(Y, contrast = "contr.sum", drop = TRUE){
  nlev <- nlevels(Y)
  lev  <- levels(Y)
  if(drop){
    X    <- model.matrix(~x,data.frame(x=factor(Y)), contrasts.arg = list(x=contrast))
    X    <- X[, -1, drop=FALSE]
  } else {
    X    <- model.matrix(~x-1,data.frame(x=factor(Y)), contrasts.arg = list(x=contrast))
  }
  attributes(X) <- list(dim = attributes(X)$dim, dimnames = attributes(X)$dimnames)
  X
}
