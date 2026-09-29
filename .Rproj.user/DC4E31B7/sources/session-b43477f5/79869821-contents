fbed.cln <- function(y, X, prior = NULL, alpha = 0.05, univ = NULL, K = 0, tol = 1e-6, maxit = 500) {

  p <- dim(X)[2]
  sig <- log(alpha)
  card <- 0
  ind <- 1:p
  sa <- pva <- sela <- NULL
  if  ( is.null( colnames(X) ) )  colnames(X) <- paste("X", 1:p, sep = "")

  runtime <- proc.time()
  if ( is.null(univ) ) {
    if ( is.null(prior) ) {
       univ <- Compositionalcln::cln.regs(y, X)
    } else  univ <- Compositionalcln::cln.condregs(y, X, prior, tol, maxit)

    stat <- univ$stat
    pval <- univ$pvalue
    n.tests <- p
  } else {
    stat <- univ$stat
    pval <- univ$pvalue
    n.tests <- 0
  }
  s <- which(pval < sig)

  if ( length(s) > 0 ) {
    sel <- which.min(pval)
    sela <- sel
    s <- s[ - which(s == sel) ]
    sa <- stat[sel]
    pva <- pval[sel]
    #########
    while ( sum(s>0) > 0 ) {
      mod <- Compositionalcln::cln.condregs(y, X[, s, drop = FALSE], prior = cbind(prior, X[, sela]), tol, maxit)
      n.tests <- n.tests + length( ind[s] )
      stat[-s] <-  -Inf
      pval[-s] <- 0
      stat[s] <- mod$stat
      pval[s] <- mod$pvalue
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel])
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0] )
      s <- s[ - which(s == sel) ]
    } ## end while ( sum(s > 0) > 0 )
    card <- sum(sela > 0)

    if ( K == 1 ) {
      mod <- Compositionalcln::cln.condregs(y, X[, -sela, drop = FALSE], prior = cbind(prior, X[, sela]), tol, maxit)
      n.tests[2] <- length( ind[-sela] )
      stat[sela] <-  -Inf
      pval[sela] <- 0
      stat[-sela] <- mod$stat
      pval[-sela] <- mod$pvalue
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel])
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      while ( sum(s>0) > 0 ) {
        mod <- Compositionalcln::cln.condregs(y, X[, s, drop = FALSE], prior = cbind(prior, X[, sela]), tol, maxit)
        n.tests[2] <- n.tests[2] + length( ind[s] )
        stat[-s] <-  -Inf
        pval[-s] <- 0
        stat[s] <- mod$stat
        pval[s] <- mod$pvalue
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel])
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
      } ## end while ( sum(s>0) > 0 )
      card <- c(card, sum(sela>0) )
    }  ## end if ( K == 1 )

    if ( K > 1 ) {
      mod <- Compositionalcln::cln.condregs(y, X[, -sela, drop = FALSE], prior = cbind(prior, X[, sela]), tol, maxit)
      n.tests[2] <- length( ind[-sela] )
      stat[sela] <-  -Inf
      pval[sela] <- 0
      stat[-sela] <- mod$stat
      pval[-sela] <- mod$pvalue
      s <- which(pval < sig)
      sel <- which.min(pval) * ( length(s)>0 )
      sa <- c(sa, stat[sel])
      pva <- c(pva, pval[sel])
      sela <- c(sela, sel[sel>0])
      s <- s[ - which(s == sel) ]
      while ( sum(s > 0) > 0 ) {
        mod <- Compositionalcln::cln.condregs(y, X[, s, drop = FALSE], prior = cbind(prior, X[, sela]), tol, maxit)
        n.tests[2] <- n.tests[2] + length( ind[s] )
        stat[-s] <-  -Inf
        pval[-s] <- 0
        stat[s] <- mod$stat
        pval[s] <- mod$pvalue
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel])
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]
      } ## end while ( sum(s>0) > 0 )
      card <- c(card, sum(sela > 0) )

      vim <- 1
      while ( vim < K  &  card[vim + 1] - card[vim] > 0  &  sela < p ) {
        vim <- vim + 1
        mod <- Compositionalcln::cln.condregs(y, X[, -sela, drop = FALSE], prior = cbind(prior, X[, sela]), tol, maxit)
        n.tests[vim + 1] <- length( ind[-sela] )
        stat[sela] <-  -Inf
        pval[sela] <- 0
        stat[-sela] <- mod$stat
        pval[-sela] <- mod$pvalue
        s <- which(pval < sig)
        sel <- which.min(pval) * ( length(s)>0 )
        sa <- c(sa, stat[sel])
        pva <- c(pva, pval[sel])
        sela <- c(sela, sel[sel>0])
        s <- s[ - which(s == sel) ]

        while ( sum(s > 0) > 0 ) {
          mod <- Compositionalcln::cln.condregs(y, X[, s, drop = FALSE], prior = cbind(prior, X[, sela]), tol, maxit)
          n.tests[vim + 1] <- n.tests[vim + 1] + length( ind[s] )
          stat[-s] <-  - Inf
          pval[-s] <- 0
          stat[s] <- mod$stat
          pval[s] <- mod$pvalue
          s <- which(pval < sig)
          sel <- which.min(pval) * ( length(s)>0 )
          sa <- c(sa, stat[sel])
          pva <- c(pva, pval[sel])
          sela <- c(sela, sel[sel>0])
          s <- s[ - which(s == sel) ]
        } ## end while ( sum(s > 0) > 0 )
        card <- c(card, sum(sela>0) )
      }  ## end while ( vim < K )
    } ## end if ( K > 1)
  } ## end if ( length(s) > 0 )

  runtime <- proc.time() - runtime
  len <- sum( sela > 0 )
  if ( len > 0 ) {
    res <- cbind(sela[1:len], sa[1:len], pva[1:len] )
    info <- matrix(nrow = length(card), ncol = 2)
    info[, 1] <- card
    info[, 2] <- n.tests
  } else {
    res <- matrix(c(0, 0, 0), ncol = 3)
    info <- matrix(c(0, p), ncol = 2)
  }
  colnames(res) <- c("Vars", "stat", "log p-value")
  rownames(info) <- paste("K=", 1:length(card)- 1, sep = "")
  colnames(info) <- c("Number of vars", "Number of tests")
  list(univ = univ, res = res, info = info, runtime = runtime)
}
