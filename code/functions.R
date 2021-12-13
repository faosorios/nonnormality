## ID: functions.R, last updated 2021-11-19, F.Osorio

## Bootstrap confidence intervals and distribution for 'eta' and 'Negentropy'
boots <- function(x, B = 1000) {
  p <- ncol(x)
  obs  <- 1:nrow(x)
  fit  <- studentFit(x[obs,], family = Student(eta = 0.25))
  eta0 <- fit$eta
  neg0 <- negentropy(eta0, p)
  stat <- matrix(0, nrow = B, ncol = 2)

  info <- fisher.matrix(fit)
  rows <- cols <- seq.int(from = p + 1, length.out = p * (p + 1) / 2)
  fisher.Scatter <- info[rows, cols]
  cols <- ncol(info)
  fisher.cross <- info[rows, cols]
  fisher.eta <- info[cols, cols]
  sigma2 <- 1 / c(fisher.eta - t(fisher.cross) %*% solve(fisher.Scatter, fisher.cross))
  c.eta <- eta0 / (1 - 2 * eta0)
  diff <- trigamma(0.5 * (1 + eta0 * p) / eta0) - trigamma(.5 / eta0)
  hdot <- 0.5 * (p * c.eta + 0.5 * (1 + eta0 * p) * diff / eta0) / eta0^2

  for (i in 1:B) {
    star <- sample(obs, replace = TRUE)
    eta  <- studentFit(x[star,], family = Student(eta = 0.25))$eta
    neg  <- negentropy(eta, p)
    stat[i,1] <- eta
    stat[i,2] <- neg
  }
  se <- sqrt(apply(stat, 2, var))

  # CI for eta
  cnames <- c("L","U")
  norm <- c(eta0 - 2 * se[1], eta0 + 2 * se[1])
  perc <- quantile(stat[,1], prob = c(0.025, 0.975))
  pivo <- 2 * eta0 - quantile(stat[,1], prob = c(0.975, 0.025))
  ci1  <- rbind(norm, perc, pivo)
  colnames(ci1) <- cnames
  rownames(ci1) <- c("normal","percentile","pivotal")

  # CI for Negentropy
  norm <- c(neg0 - 2 * se[2], neg0 + 2 * se[2])
  perc <- quantile(stat[,2], prob = c(0.025, 0.975))
  pivo <- 2 * neg0 - quantile(stat[,2], prob = c(0.975, 0.025))
  ci2  <- rbind(norm, perc, pivo)
  colnames(ci2) <- cnames
  rownames(ci2) <- c("normal","percentile","pivotal")

  o <- list(eta = eta0, neg = neg0, stat = stat, se = se, sigma2 = sigma2, var.neg = sigma2 * hdot^2, ci.eta = ci1, ci.neg = ci2)
  o
}

## Fisher information matrix for the multivariate t-distribution
fisher.matrix <- function(object) {
  ## local functions
  c.eta <- function(eta) eta / (1 - 2 * eta)
  c.phi <- function(eta, p) (1 + p * eta) / (1 + (p + 2) * eta)
  c.mu <- function(eta, p) c.phi(eta, p) / (1 - 2 * eta)
  beta.dot <- function(eta, p) {
      dif <- trigamma(.5 * (1 + p * eta) / eta) - trigamma(.5 / eta)
      -.5 * dif / eta^2
  }

  ## extract elements from 'object'
  Scatter <- object$Scatter
  eta <- object$eta

  ## checking elements
  p <- ncol(Scatter)
  if (nrow(Scatter) != p)
    stop("'Scatter' must be a square dispersion matrix")
  if (!isSymmetric(Scatter))
    Scatter <- asSymmetric(Scatter)

  ## invert and vectorizing Scatter matrix
  inv <- solve(Scatter)
  if (!isSymmetric(inv))
    inv <- asSymmetric(inv)
  vec.Scatter <- as.vector(inv)

  ## Fisher information matrix about 'center'
  fisher.center <- c.mu(eta, p) * inv
  if (!isSymmetric(fisher.center))
    fisher.center <- asSymmetric(fisher.center)

  ## Fisher information matrix about 'Scatter'
  fisher.Scatter <- 2 * c.phi(eta, p) * symm.prod(n = p, kronecker.prod(inv, inv), side = "right")
  fisher.Scatter <- fisher.Scatter + (c.phi(eta, p) - 1) * outer(vec.Scatter, vec.Scatter)
  fisher.Scatter <- .25 * dupl.cross(n = p, k = p, fisher.Scatter)
  if (!isSymmetric(fisher.Scatter))
    fisher.Scatter <- asSymmetric(fisher.Scatter)

  ## crossed Fisher information about 'Scatter' and 'eta'
  fisher.cross <- -(c.eta(eta) * (p + 2)) / ((1 + p * eta) * (1 + (p + 2) * eta))
  fisher.cross <- fisher.cross * dupl.prod(n = p, vec.Scatter, transposed = TRUE, side = "left")

  ## Fisher information about 'eta'
  dif <- trigamma(.5 / eta) - trigamma(.5 * (1 + p * eta) / eta)
  num <- 4 * (p + 2) * eta^2 - p * eta - 1
  den <- (1 - p * eta) * (1 + (p + 2) * eta)
  fisher.eta <- .25 * (dif + 2 * p * c.eta(eta) * (num / den)) / eta^4

  ## forming the Fisher information matrix
  Dp.cols <- p * (p + 1) / 2
  fisher <- matrix(0, nrow = p + Dp.cols + 1, ncol = p + Dp.cols + 1)
  fisher[1:p, 1:p] <- fisher.center
  rows <- cols <- seq.int(from = p + 1, length.out = Dp.cols)
  fisher[rows, cols] <- fisher.Scatter
  cols <- ncol(fisher)
  fisher[rows, cols] <- fisher[cols, rows] <- fisher.cross
  fisher[cols, cols] <- fisher.eta

  fisher
}

## Negentropy for a vector following a multivariate t-distribution
negentropy <- function(eta, p) {
  c.eta <- eta / (1 - 2 * eta)
  term1 <- .5 * p * (1 + log(2 * pi))
  term2 <- .5 * p * log(c.eta) - .5 * p * log(pi) + lgamma(.5 * (1 + eta * p) / eta) - lgamma(.5 / eta)
  term3 <- digamma(.5 * (1 + p * eta) / eta) - digamma(.5 / eta)
  ans <- term1 + term2 - .5 * (1 + eta * p) * term3 / eta
  ans
}

## Influence function for 'eta'
score.eta <- function(eta, p, distances) {
  c.eta <- eta / (1 - 2 * eta)
  denom <- 1 + c.eta * distances
  term1 <- log(denom)
  term2 <- digamma(.5 * (1 + p * eta) / eta) - digamma(.5 / eta)
  term3 <- (distances / denom) * (1 + p * eta) / (1 - 2 * eta) - p
  term3 <- c.eta * term3
  ans <- .5 * (term1 - term2 - term3) / eta^2
  ans
}

## Summary statistics
stats <- function(x) {
  # 1D summary stats
  xbar <- colMeans(x)
  sd <- apply(x, 2, sd)
  z <- apply(x, 2, moments) # list up forth central moments

  p <- ncol(x)
  z <- unlist(z)
  which  <- seq(p + 1, by = 5, length = p)
  skew1D <- z[which]
  which  <- seq(p + 2, by = 5, length = p)
  kurt1D <- z[which] + 3

  # multivariate summary stats
  skew <- skewness(x)
  kurt <- kurtosis(x)

  # output object
  out1D <- cbind(xbar, sd, skew1D, kurt1D)
  colnames(out1D) <- c("mean", "sd", "skewness", "kurtosis")
  z <- list(summary.1D = out1D, Scatter = cov(x), mardia = c(skew, kurt))
  z
}
