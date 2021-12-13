## ID: student.influence.R, last updated 2021-11-19, F.Osorio

influence.t <- function(x)
{ ## t-perturbation local influence
  derivatives.t <- function(x) {
    n <- dim(x)[1]
    p <- dim(x)[2]

    z <- cov.weighted(x)
    xbar <- z$mean
    S <- z$cov
    D2 <- Mahalanobis(x, xbar, S)

    ldh  <- p * (p + 3) / 2
    hess <- matrix(0, nrow = ldh, ncol = ldh)
    inv  <- solve(S)
    hess[1:p,1:p] <- n * inv
    rows <- cols <- (p + 1):ldh
    hess[rows,cols] <- .5 * n * dupl.cross(n = p, x = kronecker.prod(inv, inv))

    Delta <- matrix(0, nrow = ldh, ncol = n)
    z <- scale(x, scale = FALSE)
    Delta[1:p,] <- inv %*% t(z)
    for (i in 1:n) {
      u <- p + 2 - D2[i]
      o <- z[i,]
      Delta[1:p,i] <- u * Delta[1:p,i]
      Delta[rows,i] <- .5 * u * dupl.prod(n = p, as.vector(inv %*% outer(o, o) %*% inv - inv), transposed = TRUE)
    }
    out <- list(hess = hess, Delta = Delta)
    out
  }

  obj  <- derivatives.t(x)
  curv <- t(obj$Delta) %*% solve(obj$hess, obj$Delta)
  scaling <- sqrt(sum(diag(crossprod(curv))))
  curv <- curv / scaling

  ## compute largest eigenvectors and create the output object
  rs <- svd(curv, nv = 0)
  which <- abs(rs$d) < .Machine$double.eps
  hmax  <- rs$u[,1]
  out <- list(hmax = hmax, vectors = rs$u[,!which], total = diag(curv))

  class(out) <- "influence.t"
  out
}

print.influence.t <- function(x, idn = 3, ...)
{ # print influential observations
  hmax <- abs(x$hmax)
  cutoff <- mean(hmax) + 2 * sd(hmax)
  n <- length(hmax)
  cat("\nCutoff:", format(cutoff))

  if (is.null(idn))
    idn <- 0
  else {
    idn <- as.integer(idn)
    if(idn < 0 || idn > n)
    stop(paste("`idn' must be in { 1,..,",n,"}"))
  }
  if(idn > 0) {
    idx   <- 1:idn
    show  <- order(hmax, decreasing = TRUE)[idx]
    show  <- intersect(show, (1:n)[hmax > cutoff])
    cat("\nInfluential observations:\n")
    print(sort(show))
  }
  invisible(x)
}

plot.influence.t <- function(x, idn = 3, main = "")
{ # plotting direction of largest curvature
  hmax <- abs(x$hmax)
  cutoff <- mean(hmax) + 2 * sd(hmax)
  n <- length(hmax)

  if (is.null(idn))
    idn <- 0
  else {
    idn <- as.integer(idn)
    if(idn < 0 || idn > n)
    stop(paste("`idn' must be in { 1,..,",n,"}"))
  }
  if(idn > 0) {
    idx  <- 1:idn
    show <- order(-hmax)[idx]
    show <- intersect(show, (1:n)[hmax > cutoff])
  }

  plot(hmax, ylim = c(0,1), ylab = "hmax", main = main, font.main = 1)
  abline(h = cutoff, col = "red", lty = 2, lwd = 2)
  if (idn > 0) {
    which <- 1:n
    text(which[show], hmax[show], as.character(show), pos = 3)
  }
  invisible(x)
}
