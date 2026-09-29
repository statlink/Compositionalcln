cln.biplot <- function(y, m, S) {
  dims <- c(1, 2)
  d <- dim(S)[2]  ;   D <- d + 1
  W <- ( diag(D) - matrix(1 / D, D, D) )[, 1:d]
  e <- eigen( W %*% S %*% t(W), symmetric = TRUE)
  lam <- e$values[dims]  ;  Vp <- e$vectors[, dims, drop = FALSE]
  P <- diag(1 / sqrt(lam)) %*% t(Vp) %*% W
  G <- Vp %*% diag( sqrt(lam) )
  pve <- 100 * lam / sum( e$values[1:d] )
  labels <- colnames(y)
  if ( is.null(labels) )  labels <- paste0("Y", 1:D)

  L <- .ez(y, m, S)
  scores <- t( P %*% t(L$Z - L$M) )
  k <- max( abs(G)) / max(abs(scores) ) * 0.9
  haszero <- as.logical( Rfast::rowsums(y == 0) )

  lim <- max( abs(G), abs(scores * k) )
  plot(0, 0, type = "n", xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1,
       xlab = sprintf("PC%d (%.1f%%)", dims[1], pve[1]),
       ylab = sprintf("PC%d (%.1f%%)", dims[2], pve[2]), cex.lab = 1.3, cex.axis = 1.3)
  abline(h = 0, v = 0, lty = 3, col = "grey70")

  points(scores[!haszero, 1] * k, scores[!haszero, 2] * k, pch = 16, col = "grey40")
  points(scores[haszero, 1] * k, scores[haszero, 2] * k, pch = 16, col = "red")
  arrows(0, 0, G[, 1], G[, 2], length = 0.08, col = "steelblue")
  text(G[, 1] * 1.1, G[, 2] * 1.1, labels, col = "steelblue")
  legend("topright", pch = 16, col = c("grey40", "red"), bty = "n",
         legend = c("No zeros", "With zeros") )
}


.ez <- function(y, m, S) {
  n <- dim(y)[1]   ;   D <- dim(y)[2]  ;   d <- D - 1
  M <- matrix(m, n, d, byrow = TRUE)
  Z <- matrix(0, n, d)
  A <- rbind(diag(d), 0)
  dD <- diag(D)
  for (i in 1:n) {
    s <- which( y[i, ] > 0)
    if ( length(s) == D ) {
      Z[i, ] <- log(y[i, -D] / y[i, D])
      next
    }
    cs <- length(s) - 1
    E <- dD[s, , drop = FALSE]
    Q <- ( E[-length(s), , drop = FALSE] - matrix(E[length(s), ], cs, D, byrow = TRUE)) %*% A
    b <- log( y[i, s[-length(s)]]) - log(y[i, s[length(s)]] )
    K <- S %*% t(Q) %*% solve( Q %*% S %*% t(Q) )
    Z[i, ] <- M[i, ] + K %*% ( b - Q %*% M[i, ] )
  }
  list(Z = Z, M = M)
}


