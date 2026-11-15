cln.mle <- function(y, tol = 1e-6, maxit = 500) {
  y <- as.matrix(y)
  D <- dim(y)[2]  ;  n <- dim(y)[1]
  storage.mode(y) <- "double"
  mod <- .cln_mle_cpp(y, tol, maxit)

  n1 <- sum( Rfast::rowsums(y > 0) == D )
  y2 <- y[Rfast::rowsums(y > 0) != D, , drop = FALSE]
  mat <- matrix( as.numeric(y2 != 0), nrow = dim(y2)[1] )
  key <- apply(mat, 1, paste, collapse = ",")
  theta <- as.vector( table( factor(key, levels = unique(key)) ) )  ## same order as unique(mat)
  const <- n1 * log(n1/n) + sum( theta * log(theta/n) )
  patterns <- cbind(unique(mat), theta/n)

  m <- as.vector(mod$m)
  mesi <- c(exp(m), 1)  ;  mesi <- mesi/sum(mesi)
  list( m = m, mesi = mesi, S = mod$S, loglik = mod$loglik + const, patterns = patterns,
        iters = mod$iters )
}
