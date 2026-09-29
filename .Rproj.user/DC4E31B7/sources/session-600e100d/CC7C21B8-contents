cln.pca <- function(y, m, S) {
  d <- dim(S)[2]  ;   D <- d + 1
  W <- ( diag(D) - matrix(1 / D, D, D) )[, 1:d]
  e <- eigen( W %*% S %*% t(W), symmetric = TRUE )
  lam <- e$values[1:d]
  V <- e$vectors[, 1:d, drop = FALSE]
  labels <- colnames(y)
  if ( is.null(labels) )  labels <- paste0("Y", 1:D)
  rownames(V) <- labels  ;  colnames(V) <- paste0("PC", 1:d)
  L <- .ez(y, m, S)
  scores <- ( L$Z - L$M ) %*% t(W) %*% V
  colnames(scores) <- paste0("PC", 1:d)
  list( values = lam, prop = lam / sum(lam), loadings = V, scores = scores )
}
