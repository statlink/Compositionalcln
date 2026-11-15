## Log-likelihood of the CLN regression model, used by cln.regs() and cln.condregs().
.loglik.zero.norm.reg <- function(B, S, x1, x2, full, y1, y2, sly1, Q.list, obs.list) {
  n1 <- dim(y1)[1]  ;  n2 <- dim(y2)[1]
  Sinv <- solve(S)
  res1 <- full - x1 %*% B
  quad1 <- sum( Rfast::rowsums( (res1 %*% Sinv) * res1 ) )
  ll1 <-  -0.5 * n1 * log( det(2 * pi * S) ) - 0.5 * quad1 - sly1

  ll2 <- 0
  for ( i in seq_len(n2) ) {
    Qi <- Q.list[[ i ]]
    b <- obs.list[[ i ]]
    mu_i <- drop( crossprod(B, x2[i, ]) )
    muA <- drop( Qi %*% mu_i )
    SA  <- Qi %*% S %*% t(Qi)
    Ci <- length(b) + 1
    z1 <- y2[i, ]  ;  z1 <- z1[z1 > 0]
    if ( Ci > 2 ) {
      ll2 <- ll2 - 0.5 * log( det(2 * pi * SA) ) - 0.5 * t(b - muA) %*% solve(SA, b - muA) - sum( log(z1) )
    } else {
      ll2 <- ll2 - 0.5 * log(2 * pi * SA) - 0.5 * (b - muA)^2 / SA - sum( log(z1) )
    }
  }
  as.numeric(ll1 + ll2)
}
