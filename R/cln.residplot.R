cln.residplot <- function(y, m, S, x = NULL, beta = NULL) {

  n <- dim(y)[1]   ;   D <- dim(y)[2]   ;   d <- D - 1
  if ( is.null(x) ) {
    M <- matrix(m, n, d, byrow = TRUE)
    p <- 1
  } else {
    x <- model.matrix(y~., data = as.data.frame(x) )
    M <- x %*% beta
    p <- dim(x)[2]
  }
  A <- rbind(diag(d), 0)  ;  dD <- diag(D)
  m2 <- rep(NA, n)
  dof <- Rfast::rowsums(y > 0) - 1
  haszero <- dof < d

  for ( i in which(dof >= 1) ) {
    s <- which(y[i, ] > 0)  ;  cs <- length(s) - 1
    E <- dD[s, , drop = FALSE]
    Q <- ( E[-length(s), , drop = FALSE] - matrix(E[length(s), ], cs, D, byrow = TRUE) ) %*% A
    b <- log( y[i, s[-length(s)]] ) - log( y[i, s[length(s)]] )
    res <- b - Q %*% M[i, ]
    m2[i] <- drop( crossprod(res, solve(Q %*% S %*% t(Q), res)) )
  }
  lp <- pchisq(m2, dof, log.p = TRUE)
  z <- qnorm(lp, log.p = TRUE)

  qq <- qqnorm(z, plot.it = FALSE)
  plot(qq, pch = 16, col = ifelse(haszero, "red", "grey40"), cex.axis = 1.3, cex.lab = 1.3,
       xlab = "Theoretical N(0,1) quantiles", ylab = "Normal scores of the p-values")
  abline(0, 1)
  legend("topleft", pch = 16, col = c("grey40", "red"), bty = "n",
         legend = c("No zeros", "With zeros"))
  list(m2 = m2, dof = dof, cdf = exp(lp), z = z)
}

