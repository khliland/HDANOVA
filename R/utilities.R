# PCA for internal use
.pca <- function(X, ncomp, scale=FALSE, ...){
  X <- as.matrix(unclass(X))
  if(!inherits(X,'matrix'))
    stop("'X' must be a matrix")
  if(missing(ncomp))
    ncomp <- min(c(nrow(X)-1, ncol(X)))
  y   <- structure(matrix(rnorm(nrow(X)), ncol=1), rownames = rownames(X))
  dat <- data.frame(y=y, X = I(X))
  PCR <- pls::pcr(y ~ X, ncomp = ncomp, data = dat, scale = scale, ...)
  mod <- list(loadings=PCR$loadings, scores=PCR$scores, Xmeans=PCR$Xmeans, explvar=PCR$Xvar/PCR$Xtotvar*100, PCA = PCR)
  attr(mod$loadings, 'explvar') <- attr(mod$scores, 'explvar') <- mod$explvar
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

