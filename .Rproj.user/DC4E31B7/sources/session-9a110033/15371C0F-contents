cln.regs <- function(y, X, tol = 1e-6, maxit = 500) {

  D <- dim(y)[2]  ;  d <- D - 1  ;  n <- dim(y)[1]
  p <- dim(X)[2]
  stat <- pvalue <- numeric(p)

  sel1 <- Rfast::rowsums(y > 0) == D
  y1 <- y[sel1, ]                       ;  n1 <- dim(y1)[1]
  X1 <- X[sel1, , drop = FALSE]
  sly1 <- sum( log(y1) )

  y2 <- y[!sel1, ]
  X2 <- X[!sel1, , drop = FALSE]
  if (n - n1 == 1)  y2 <- matrix(y2, nrow = 1)
  n2 <- dim(y2)[1]
  full <- log( y1[, -D] / y1[, D] )
  mat <- matrix( as.numeric(y2 != 0), nrow = n2 )

  theta <- table( apply(1 - mat, 1, paste, collapse = ",") )
  theta <- as.vector(theta)
  const <- n1 * log(n1/n) + sum( theta * log(theta/n) )
  patterns <- cbind(unique(mat), theta/n)

  F <- function(d)  cbind( diag(d), -1)
  H <- function(d)  diag(d) + 1
  com <- t( F(d) ) %*% solve( H(d) )

  Q.list <- vector("list", n2)  ;  obs.list <- vector("list", n2)
  for ( i in 1:n2 ) {
    z <- mat[i, ]  ;  C <- sum(z)  ;  c <- C - 1
    Sm <- diag(z)
    Sm <- matrix( Sm[Rfast::rowsums(Sm) > 0, ], ncol = D)
    Q.list[[ i ]] <- F(c) %*% Sm %*% com
    z1 <- y2[i, ]  ;  z1 <- matrix( z1[ z1 > 0 ], nrow = 1 )
    obs.list[[ i ]] <- if ( C > 2 )  drop( log(z1[, -C] / z1[, C]) )  else log(z1[, 1] / z1[, 2])
  }

  lik0 <- Compositionalcln::cln.mle(y, tol, maxit)$loglik 

  for ( j in 1:p ) {
    x1 <- cbind(1, X1[, j])
    x2 <- cbind(1, X2[, j])
    B <- mziln::mziln(y, X[, j])$be
    res1 <- full - x1 %*% B
    S <- crossprod(res1) / n1
    Ez <- matrix(0, n2, d)
    loglik.old <- .loglik.zero.norm.reg(B, S, x1, x2, full, y1, y2, sly1, Q.list, obs.list)
    for ( it in 1:maxit ) {

      EzzVarSum <- 0
      for ( i in 1:n2 ) {
        Qi <- Q.list[[ i ]]
        b  <- obs.list[[ i ]]
        mu_i <- drop( crossprod(B, x2[i, ]) )
        muA <- drop( Qi %*% mu_i )
        SA  <- Qi %*% S %*% t(Qi)
        K   <- S %*% t(Qi) %*% solve(SA)
        ez  <- mu_i + drop( K %*% (b - muA) )
        vz  <- S - K %*% Qi %*% S
        Ez[i, ] <- ez
        EzzVarSum <- EzzVarSum + vz
      }
      Z <- rbind(full, Ez)
      XX <- rbind(x1, x2)
      B <- solve( crossprod(XX), crossprod(XX, Z) )         # GLS reduces to OLS (common S)
      res1 <- full - x1 %*% B
      res2 <- Ez - x2 %*% B
      S <- ( crossprod(res1) + crossprod(res2) + EzzVarSum ) / n
      loglik.new <- .loglik.zero.norm.reg(B, S, x1, x2, full, y1, y2, sly1, Q.list, obs.list)
      if ( abs(loglik.new - loglik.old) < tol )  break
      loglik.old <- loglik.new
    }
    lik1 <- loglik.new + const
    stat[j] <- 2 * (lik1 - lik0)
    dm <- dim(B)
    dof <- (dm[1] - 1) * dm[2]
    pvalue[j] <- pchisq(stat[j], dof, lower.tail = FALSE, log.p = TRUE)
  }  ## end  for ( j in 1:p ) {

  list(stat = stat, pvalue = pvalue)
}
