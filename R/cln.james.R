cln.james <- function(y1, y2, a = 0.05, tol = 1e-6, maxit = 500) {

  d <- dim(y1)[2] - 1  ## dimensionality of the data
  n1 <- dim(y1)[1]   ;   n2 <- dim(y2)[1]  ## sample sizes
  mod1 <- Compositionalcln::cln.mle(y1, tol, maxit)
  Compositionalcln::mod2 <- cln.mle(y2, tol, maxit)
  ybar1 <- mod1$m   ;   ybar2 <- mod2$m
  dbar <- ybar2 - ybar1  ## difference of the two mean vectors
  A1 <- mod1$S/n1   ;   A2 <- mod2$S/n2
  V <- A1 + A2  ## covariance matrix of the difference
  Vinv <- solve(V)
  test <- sum( dbar %*% Vinv * dbar )
  b1 <- Vinv %*% A1   ;   b2 <- Vinv %*% A2
  trb1 <- sum( diag(b1) )    ;   trb2 <- sum( diag(b2) )

  A <- 1 + ( trb1^2/(n1 - 1) + trb2^2/(n2 - 1) ) / (2 * p)
  B <- ( sum(b1^2) / (n1 - 1) + sum(b2^2)/(n2 - 1) + 0.5 * trb1 ^ 2/ (n1 - 1) + 0.5 * trb2^2/(n2 - 1) ) / (p * (p + 2))
  x2 <- qchisq(1 - a, d)
  delta <- (A + B * x2)
  twoha <- x2 * delta  ## corrected critical value of the chi-square
  pvalue <- pchisq(test/delta, d, lower.tail = FALSE)  ## p-value of the test statistic
  info <- c(test, pvalue, delta, twoha)
  names(info) <- c("test", "p-value", "correction", "corrected.critical")
  info
}
