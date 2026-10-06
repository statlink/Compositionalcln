boot.clnmle <- function(y, tol = 1e-6, maxit = 500, R = 1000) {
  D <- dim(y)[2]  ;  d <- D - 1
  n <- dim(y)[1]
  y1 <- y[Rfast::rowsums(y > 0) == D, ] ;  n1 <- dim(y1)[1]
  sly1 <- sum( log(y1) )
  y2 <- y[Rfast::rowsums(y > 0) != D, ]
  if (n - n1 == 1)  y2 <- matrix(y2, nrow = 1)  ;  n2 <- dim(y2)[1]
  sly2 <- sum( log(y2[y2>0]) )

  full <- log( y1[, -D] / y1[, D] )
  mfull <- Rfast::colsums(full)
  sfull <- crossprod(full)
  mat <- matrix( as.numeric(y2 != 0), nrow = n2 )

  F <- function(d)  cbind( diag(d), -1)
  H <- function(d)  diag(d) + 1
  m <- mfull/n1   ;   S <- ( (n1 - 1)/n1 ) * var(full)
  Q.list <- vector("list", n2)  ;  obs.list <- vector("list", n2)
  com <- t( F(d) ) %*% solve( H(d) )

  for ( i in 1:n2 ) {
    z <- mat[i, ]  ;  C <- sum(z)  ;  c <- C - 1
    Sm <- diag(z)
    Sm <- matrix( Sm[Rfast::rowsums(Sm) > 0, ], ncol = D )
    Q.list[[ i ]] <- F(c) %*% Sm %*% com
    z1 <- y2[i, ]  ;  z1 <- z1[ z1 > 0 ]
    obs.list[[ i ]] <- if ( C > 2 )  log(z1[-C] / z1[C])  else log(z1[1] / z1[2])
  }

  Ez <- matrix(0, n2, d)
  bm <- matrix(nrow = R, ncol = d)
  bmesi <- matrix(nrow = R, ncol = D)
  bS <- list()

  S <- matrix(NA, nrow = d, ncol = d)
  bS <- as.list( rep(list( S), R) )

  for ( vim in 1:R ) {
    id1 <- rangen::Sample.int(n1, n1, replace = TRUE)
    yb1 <- y1[id1, ]
    fullb <- full[id1, ]
    id2 <- rangen::Sample.int(n2, n2, replace = TRUE)
    yb2 <- y2[id2, ]
    fullb <- full[id1, ]
    slyb1 <- sum( log(yb1) )
    slyb2 <- sum( log(yb2[yb2>0]) )
    mfullb <- Rfast::colsums(fullb)
    sfullb <- crossprod(fullb)
    m <- mfullb/n1   ;   S <- ( (n1 - 1)/n1 ) * var(fullb)
    Q.listb <- Q.list[id2]
    obs.listb <- obs.list[id2]
    loglik.old <- .loglik.zero.norm(m, S, full, yb1, yb2, slyb1, slyb2, Q.listb, obs.listb)
    for ( it in 1:maxit ) {
      EzzSum <- 0
      for ( i in 1:n2 ) {
        Qi <- Q.listb[[ i ]]
        b <- obs.listb[[ i ]]
        muA <- drop(Qi %*% m)
        SA <- Qi %*% S %*% t(Qi)
        K <- S %*% t(Qi) %*% solve(SA)
        ez <- drop( m + K %*% (b - muA) )
        vz <- S - K %*% Qi %*% S
        Ez[i, ] <- ez
        EzzSum <- EzzSum + vz + tcrossprod(ez)
      }
      m <- ( mfullb + Rfast::colsums(Ez) ) / n
      S <- ( sfullb + EzzSum ) / n - tcrossprod(m)
      loglik.new <- .loglik.zero.norm(m, S, full, y1, y2, sly1, sly2, Q.list, obs.list)
      if ( abs(loglik.new - loglik.old) < tol )  break
      loglik.old <- loglik.new
    }
    bm[vim, ] <- m
    mesi <- c(exp(m), 1)  ;  bmesi[vim, ] <- mesi/sum(mesi)
    bS[[ vim ]] <- S
  }

  list(m = bm, mesi = bmesi, S = bS )
}


