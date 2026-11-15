boot.clnreg <- function(y, x, tol = 1e-6, maxit = 500, R = 1000) {
  y <- as.matrix(y)
  D <- dim(y)[2]  ;  d <- D - 1
  x <- model.matrix( y~., data = as.data.frame(x) )
  storage.mode(y) <- "double"
  mod <- .cln_bootreg_cpp(y, x, tol, maxit, R)
  if ( mod$fails > 0 )
    warning( mod$fails, " of the ", R, " bootstrap fits failed and are excluded from sigma." )

  bbe <- mod$beta
  sigma <- cov(bbe, use = "complete.obs")
  namx <- colnames(x)
  namy <- colnames(y)
  if ( is.null(namy) ) {
    namy <- paste("Y", 1:d, sep = "")
  } else namy <- namy[-D]
  ## the rows of bbe are vec(B): the covariates vary fastest within each log-ratio
  nam <- paste( rep(namy, each = length(namx)), namx, sep = ":" )
  colnames(sigma) <- rownames(sigma) <- colnames(bbe) <- nam
  list(beta = bbe, sigma = sigma)
}
