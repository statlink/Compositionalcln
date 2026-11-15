cln.condregs <- function(y, X, prior, tol = 1e-6, maxit = 500) {
  y <- as.matrix(y)  ;  X <- as.matrix(X)
  prior <- model.matrix( ~., data = data.frame(prior) )
  prior <- as.matrix(prior[, -1, drop = FALSE])
  storage.mode(y) <- "double"  ;  storage.mode(X) <- "double"  ;  storage.mode(prior) <- "double"
  mod <- .cln_condregs_cpp(y, X, prior, tol, maxit)
  if ( mod$fails > 0 )
    warning( mod$fails, " of the ", dim(X)[2], " regressions failed (e.g. a constant or ",
             "collinear predictor) and are returned as NA." )
  list( stat = as.vector(mod$stat), pvalue = as.vector(mod$pvalue) )
}
