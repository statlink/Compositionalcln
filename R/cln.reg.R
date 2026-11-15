cln.reg <- function(y, x, tol = 1e-6, maxit = 500, seb = FALSE, xnew = NULL) {
  y <- as.matrix(y)
  D <- dim(y)[2]  ;  d <- D - 1  ;  n <- dim(y)[1]
  x <- model.matrix( y~., data = as.data.frame(x) )
  p <- dim(x)[2]
  storage.mode(y) <- "double"
  mod <- .cln_reg_cpp(y, x, tol, maxit, seb)

  sel1 <- Rfast::rowsums(y > 0) == D
  n1 <- sum(sel1)
  mat <- matrix( as.numeric(y[!sel1, , drop = FALSE] != 0), ncol = D )
  key <- apply(mat, 1, paste, collapse = ",")
  theta <- as.vector( table( factor(key, levels = unique(key)) ) )  ## same order as unique(mat)
  const <- n1 * log(n1/n) + sum( theta * log(theta/n) )
  patterns <- cbind(unique(mat), theta/n)

  B <- mod$beta
  colnames(B) <- colnames(y[, -D, drop = FALSE])
  rownames(B) <- colnames(x)

  est <- NULL
  if ( !is.null(xnew) ) {
    xnew <- model.matrix(~., as.data.frame(xnew) )
    ma <- cbind( exp( xnew %*% B ), 1 )
    est <- ma / Rfast::rowsums(ma)  ## fitted values
    colnames(est) <- colnames(y)
  }

  loglik <- mod$loglik + const
  BIC <- -2 * loglik + ( p * d + 0.5 * d * (d + 1) + dim(patterns)[1] ) * log(n)
  out <- list(patterns = patterns, iters = mod$iters, loglik = loglik, BIC = BIC, S = mod$S,
              beta = B, est = est)
  if ( seb ) {
    se.beta <- mod$se.beta
    dimnames(se.beta) <- dimnames(B)
    out$se.beta <- se.beta
  }
  out
}
