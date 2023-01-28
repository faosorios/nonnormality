# Mardia multivariate skewness
mardia1D <- function(x, y) {
  # x: skewness parameter
  # y: degrees of freedom
  skew <- x
  df <- y
  a <- sqrt(df / pi) * (gamma((df - 1) / 2) / gamma(df / 2))
  delta <- skew / sqrt(1 + skew^2)
  val <- a * df * delta * (3 - delta^2) / (df - 3)
  val^2
}

mardia <- function(x, y, p = 5) {
  # x: skewness parameter
  # y: degrees of freedom
  skew <- x
  df <- y
  ones <- rep(1, p)
  Omega <- diag(p)
  a <- sqrt(df / pi) * (gamma((df - 1) / 2) / gamma(df / 2))
  delta <- (skew / sqrt(1 + p * skew^2)) * ones
  s <- kronecker.prod(delta, Omega) + outer(vec(Omega), delta) + kronecker.prod(diag(p), delta) - kronecker.prod(delta, outer(delta, delta))
  s <- (a * df / (df - 3)) * s
  val <- sum(diag(crossprod(s)))
  val
}

df <- seq(from = 3.01, to = 3.1, length = 30)
skew <- seq(from = -2, to = 2, length = 30)
z <- outer(skew, df, mardia1D)
z5 <- matrix(0, nrow = 30, ncol = 30)
for (i in 1:30) {
  for (j in 1:30)
    z5[i,j] <- mardia(skew[i], df[j], p = 5)
}
z25 <- matrix(0, nrow = 30, ncol = 30)
for (i in 1:30) {
  for (j in 1:30)
    z25[i,j] <- mardia(skew[i], df[j], p = 25)
}
