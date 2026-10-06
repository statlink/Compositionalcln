cln.bcr <- function(y, tol = 1e-6, maxit = 500, alpha = 0.05, R = 1000, dg = TRUE, hg = TRUE, ploty = TRUE) {

  mesi <- Compositionalcln::cln.mle(y, tol = tol, maxit = maxit)$mesi
  bmesi <- Compositionalcln::boot.clnmle(y, tol = tol, maxit = maxit, R = R)$mesi

  fit <- MASS::cov.rob(bmesi[, 1:2], quantile.used = ceiling( (1 - alpha) * R), method = "mve")
  cr <- predict( cluster::ellipsoidhull( bmesi[fit$best, 1:2] ) )
  cr[cr < 0] <- 0
  cr[cr > 1] <- 1
  cr3 <- pmax(0, 1 - Rfast::rowsums(cr) )
  cr <- cbind(cr, cr3)
  confr <- cr / Rfast::rowsums(cr)

  nam <- paste("Y", 1:3, sep = "")

  b1 <- c(0.5, 0, 1, 0.5)
  b2 <- c(sqrt(3)/2, 0, 0, sqrt(3)/2)
  b <- cbind(b1, b2)
  plot(b[, 1], b[, 2], type = "l", xlab = " ", ylab = " ", pty = "s",
       xaxt = "n", yaxt = "n", bty = "n", lwd = 2)
  proj <- matrix(c(0, 1, 0.5, 0, 0, sqrt(3)/2), ncol = 2)
  Y <- rbind(y, mesi) %*% proj
  ellipse <- confr %*% proj
  n <- dim(y)[1]
  symb <- rep(16, n)

  ind <- which(y == 0, arr.ind = TRUE)[, 1]
  if ( length(ind) > 0 ) {
    symb[ind] <- 4
  }

  text( b[1, 1], b[1, 2] + 0.02, nam[3], col = "black", font = 2 )
  text( b[2, 1] + 0.02, b[2, 2] - 0.02, nam[1], col = "black", font = 2 )
  text( b[3, 1] - 0.02, b[2, 2] - 0.02, nam[2], col = "black", font = 2 )

  if ( dg ) {
    a1 <- matrix(0, nrow = 11, ncol = 3)
    a1[, 2] <- seq(0, 1, by = 0.1)
    a1[, 3] <- seq(1, 0, by = -0.1)
    ## right
    b1 <- a1 %*% proj
    ## left
    a2 <- cbind(a1[, 2], a1[, 1], a1[, 3])
    b2 <- a2 %*% proj
    ## horizontal
    a3 <- cbind(a1[, 2], a1[, 3], a1[, 1])
    b3 <- a3 %*% proj

    for ( i in 2:dim(b1)[1] ) {
      segments(x0 = b1[i, 1], y0 = b1[i, 2], x1 = b3[12 - i, 1], y1 = b3[i, 2], col = "lightgrey", lty = 2)
    }
    for (i in 1:(dim(b1)[1] - 1 ) ) {
      segments(x0 = b2[i, 1], y0 = b2[i, 2], x1 = b3[i, 1], y1 = b3[12 - i, 2], col = "lightgrey", lty = 2)
    }
    lines(b[, 1], b[, 2])
    for ( i in 2:( dim(b1)[1] - 1 ) )  {
      text(b1[i, 1] + 0.025, b1[i, 2] + 0.025, a1[i, 3], cex = 1)
      text(b2[i, 1] - 0.025, b2[i, 2] + 0.02, a1[i, 2], cex = 1)
      text(b3[i, 1], b3[i, 2] - 0.02, a1[i, 3], cex = 1)
    }
  } ## end if dg

  if ( hg ) {
    a1 <- matrix(0, nrow = 11, ncol = 3)
    a1[, 2] <- seq(0, 1, by = 0.1)
    a1[, 3] <- seq(1, 0, by = -0.1)
    ## right
    b1 <- a1 %*% proj
    ## left
    a2 <- cbind(a1[, 2], a1[, 1], a1[, 3])
    b2 <- a2 %*% proj
    for ( i in 2:c( dim(b1)[1] - 1) ) {
      segments(x0 = b1[i, 1], y0 = b1[i, 2], x1 = b2[i, 1], y1 = b2[i, 2], col = "lightgrey", lty = 2)
    }
    lines(b[, 1], b[, 2])
    for ( i in 2:( dim(b1)[1] - 1 ) )  {
      text(b1[i, 1] + 0.025, b1[i, 2] + 0.025, a1[i, 3], cex = 1)
      text(b2[i, 1] - 0.025, b1[i, 2] + 0.02, a1[i, 2], cex = 1)
    }
  } ## end if hg

  if ( ploty )  points( Y[1:n, 1], Y[1:n, 2], pch = symb, lwd = 2 )
  points( Y[c(n + 1), 1], Y[c(n + 1), 2], pch = 20, col = 4 )
  lines(ellipse, col = 4)

}

