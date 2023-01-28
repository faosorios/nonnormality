## ID: boots.R, last updated 2021-11-19, F.Osorio

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

## Bootstrap confidence intervals, skew-normal model
boots.sn <- function(x, B = 1000, nsim = 10000) {
  negentropy <- function(Sigma, Omega, alpha, nsim) {
    term1 <- c(determinant(Sigma, log = TRUE)$modulus)
    term2 <- c(determinant(Omega, log = TRUE)$modulus)
    tau <- sqrt(c(crossprod(alpha, Omega %*% alpha)))
    w <- as.vector(rsn(nsim, xi = 0, omega = 1, alpha = tau))
    ans <- .5 * term1 - .5 * term2 + mean(log(2 * pnorm(tau * w)))
    ans
  }

  p <- ncol(x)
  obs <- 1:nrow(x)
  fn0 <- studentFit(x[obs,], family = Student(eta = 0))
  fsn <- msnFit(x[obs,], description = "skew normal")
  Sigma0 <- fn0$Scatter
  Omega0 <- fsn@fit$estimated$Omega
  alpha0 <- fsn@fit$estimated$alpha
  neg0 <- negentropy(Sigma0, Omega0, alpha0, nsim)
  stat <- rep(0, B)

  cat(" Progress:\n")
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (i in 1:B) {
    star <- sample(obs, replace = TRUE)
    Sigma <- studentFit(x[star,], family = Student(eta = 0))$Scatter
    fm <- msnFit(x[star,], description = "skew normal")
    Omega <- fm@fit$estimated$Omega
    alpha <- fm@fit$estimated$alpha
    neg  <- negentropy(Sigma, Omega, alpha, nsim)
    stat[i] <- neg
    setTxtProgressBar(pb, i)
  }
  se <- sqrt(var(stat))
  close(pb)

  # CI for Negentropy
  cnames <- c("L","U")
  norm <- c(neg0 - 2 * se, neg0 + 2 * se)
  perc <- quantile(stat, prob = c(0.025, 0.975))
  pivo <- 2 * neg0 - quantile(stat, prob = c(0.975, 0.025))
  ci <- rbind(norm, perc, pivo)
  colnames(ci) <- cnames
  rownames(ci) <- c("normal","percentile","pivotal")

  o <- list(neg = neg0, stat = stat, se = se, CI = ci)
  o
}

## Bootstrap confidence intervals, skew-t model
boots.st <- function(x, B = 1000, nsim = 10000) {
  negentropy <- function(Sigma, Omega, alpha, nu, nsim) {
    p <- ncol(Sigma)
    term1 <- c(determinant(Sigma, log = TRUE)$modulus) + p * (1 + log(2 * pi))
    term2 <- c(determinant(Omega, log = TRUE)$modulus)
    term3 <- lgamma(.5 * (nu + p)) - lgamma(.5 * nu) - .5 * p * log(nu * pi)
    term4 <- .5 * (nu + p) * (digamma(.5 * (nu + p)) - digamma(.5 * nu))
    tau <- sqrt(c(crossprod(alpha, Omega %*% alpha)))
    w <- as.vector(rst(nsim, xi = 0, omega = 1, alpha = tau, nu = nu + p - 1))
    u <- sqrt(nu + p) * w / sqrt(nu + p - 1 + w^2)
    ans <- .5 * term1 - .5 * term2 + term3 - term4 + mean(log(2 * pt(tau * u, df = nu + p)))
    ans
  }

  p <- ncol(x)
  obs <- 1:nrow(x)
  fn0 <- studentFit(x[obs,], family = Student(eta = 0))
  fst <- mstFit(x[obs,], description = "skew normal")
  Sigma0 <- fn0$Scatter
  Omega0 <- fst@fit$estimated$Omega
  alpha0 <- fst@fit$estimated$alpha
  nu0 <- fst@fit$estimated$nu
  neg0 <- negentropy(Sigma0, Omega0, alpha0, nu0, nsim)
  stat <- matrix(0, nrow = B, ncol = 2)

  cat(" Progress:\n")
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for (i in 1:B) {
    star <- sample(obs, replace = TRUE)
    Sigma <- studentFit(x[star,], family = Student(eta = 0))$Scatter
    fm <- mstFit(x[star,], description = "skew-t")
    Omega <- fm@fit$estimated$Omega
    alpha <- fm@fit$estimated$alpha
    nu <- fm@fit$estimated$nu
    neg  <- negentropy(Sigma, Omega, alpha, nu, nsim)
    stat[i,1] <- nu
    stat[i,2] <- neg
    setTxtProgressBar(pb, i)
  }
  se <- rep(0, 2)
  se[1] <- sqrt(var(stat[,1], na.rm = TRUE))
  se[2] <- sqrt(var(stat[,2], na.rm = TRUE))
  #se <- sqrt(apply(stat, 2, var))
  close(pb)

  # CI for eta
  cnames <- c("L","U")
  norm <- c(nu0 - 2 * se[1], nu0 + 2 * se[1])
  perc <- quantile(stat[,1], prob = c(0.025, 0.975), na.rm = TRUE)
  pivo <- 2 * nu0 - quantile(stat[,1], prob = c(0.975, 0.025), na.rm = TRUE)
  ci1  <- rbind(norm, perc, pivo)
  colnames(ci1) <- cnames
  rownames(ci1) <- c("normal","percentile","pivotal")

  # CI for Negentropy
  norm <- c(neg0 - 2 * se[2], neg0 + 2 * se[2])
  perc <- quantile(stat[,2], prob = c(0.025, 0.975), na.rm = TRUE)
  pivo <- 2 * neg0 - quantile(stat[,2], prob = c(0.975, 0.025), na.rm = TRUE)
  ci2  <- rbind(norm, perc, pivo)
  colnames(ci2) <- cnames
  rownames(ci2) <- c("normal","percentile","pivotal")

  o <- list(nu = nu0, neg = neg0, stat = stat, se = se, CI.nu = ci1, CI.neg = ci2)
  o
}

