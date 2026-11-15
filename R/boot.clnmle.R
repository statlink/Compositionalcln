boot.clnmle <- function(y, tol = 1e-6, maxit = 500, R = 1000) {
  y <- as.matrix(y)
  D <- dim(y)[2]
  storage.mode(y) <- "double"
  mod <- .cln_bootmle_cpp(y, tol, maxit, R)
  if ( mod$fails > 0 )
    warning( mod$fails, " of the ", R, " bootstrap fits failed and are returned as NA." )

  m <- mod$m
  ma <- cbind( exp(m), 1 )
  mesi <- ma / Rfast::rowsums(ma)
  bS <- lapply( seq_len(R), function(i) mod$S[, , i] )
  nam <- colnames(y)
  if ( !is.null(nam) ) {
    colnames(m) <- nam[-D]  ;  colnames(mesi) <- nam
    bS <- lapply( bS, function(s) { dimnames(s) <- list(nam[-D], nam[-D])  ;  s } )
  }
  list(m = m, mesi = mesi, S = bS)
}
