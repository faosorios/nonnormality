## ID: stats.R, last updated 2021-11-19, F.Osorio

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
