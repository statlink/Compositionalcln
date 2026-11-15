cln.regs <- function(y, X, tol = 1e-6, maxit = 500) {
  y <- as.matrix(y)  ;  X <- as.matrix(X)
  storage.mode(y) <- "double"  ;  storage.mode(X) <- "double"
  mod <- .cln_regs_cpp(y, X, tol, maxit)
  if ( mod$fails > 0 )
    warning( mod$fails, " of the ", dim(X)[2], " regressions failed (e.g. a constant ",
             "predictor) and are returned as NA." )
  list( stat = as.vector(mod$stat), pvalue = as.vector(mod$pvalue) )
}
