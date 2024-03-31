gpca <- function(xcs, states, pcs = 1:min(5,ncol(xcs)), tol = 1e-15) {
  # Group-wise Principal Component Analysis.
  #
  # Args:
  # xcs: [NxM] preprocessed billinear data set
  # states: {S} List with the groups of variables.
  # pcs: [1xA] Principal Components considered (e.g., pcs = 1:2 selects the
  #     first two PCs). By default, pcs = 1:ncol(xcs)
  # tol: [1x1] tolerance value
  #
  # Returns:
  # p: [MxA] matrix of loadings.
  # t: [NxA] matrix of scores.
  # bel: [Ax1] correspondence between PCs and States.
  # e: [NxM] matrix of residuals.

  # Set default values

  # Preprocessing
  pcs <- unique(pcs)
  pcs <- pcs[pcs != 0]
  #pcs <- pcs[pcs < qr.R(qr(xcs))$rank]
  A <- length(pcs)

  # Validate values of input data
  stopifnot(is.list(states))
  for (i in seq_along(states)) {
    stopifnot(all(states[[i]] > 0) && all(states[[i]] %% 1 == 0))
    stopifnot(all(states[[i]] <= ncol(xcs)))
  }
  stopifnot(all(pcs > 0))

  # Main code
  map <- crossprod(xcs) #t(xcs) %*% xcs
  I <- diag(nrow(map))
  B <- I

  p <- matrix(NA, ncol(xcs), A)
  t <- matrix(NA, nrow(xcs), A)
  bel <- integer(A)
  for (j in seq_len(max(pcs))) {
    R <- matrix(0, ncol(xcs), length(states))
    S <- matrix(0, nrow(xcs), length(states))
    for (i in seq_along(states)) {
      map_aux <- matrix(0, nrow(map), ncol(map))
      map_aux[states[[i]], states[[i]]] <- map[states[[i]], states[[i]]]
      if (any(map_aux > tol)) {
        ev <- eigen(map_aux)
        ind <- which.max(diag(ev$values))
        V <- ev$vectors
#        R[, i] <- V[, ind, drop=FALSE] / sqrt(t(V[, ind, drop=FALSE]) %*% B %*% V[, ind, drop=FALSE])
        R[, i] <- V[, ind, drop=FALSE] %*% pracma::pinv(sqrt(t(V[, ind, drop=FALSE]) %*% B %*% V[, ind, drop=FALSE]))
        S[, i] <- xcs %*% R[, i, drop=FALSE]
      }
    }

    sS <- colSums(S ^ 2)
    ind <- which.max(sS)

    q <- B %*% R[, ind, drop=FALSE]
    p[, j] <- R[, ind, drop=FALSE] / norm(R[, ind, drop=FALSE])
    t[, j] <- xcs %*% p[, j]
    bel[j] <- ind
    Iq <- I - tcrossprod(q)
    map <- Iq %*% map %*% Iq # %*% t(q))
    xcs <- xcs %*% Iq #q %*% t(q))
    B <- B %*% Iq # %*% t(q))
    #    map <- (I - q %*% t(q)) %*% map %*% (I - tcrossprod(q)) # %*% t(q))
    #    xcs <- xcs %*% (I - tcrossprod(q)) #q %*% t(q))
    #    B <- B %*% (I - tcrossprod(q)) # %*% t(q))
  }

  # Postprocessing
  p <- p[, pcs]
  t <- t[, pcs]
  bel <- bel[pcs]

  e <- xcs

  return(list(loadings = p, scores = t, selBlock = bel, residual = e))
}

