## ID: envelope.R, last updated 2022-08-22, F.Osorio

envel.sn <- function(object, reps = 50, conf = 0.95, plot.it = TRUE, main = "Transformed distances Q-Q plot")
{ ## simulated envelope
  envel <- function(n, mu, Omega, alpha, reps, conf) {
    conf <- 1 - conf
    p <- length(mu)
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    elims <- matrix(0, nrow = n, ncol = reps)
    for (i in 1:reps) {
      x <- rmvsnorm(n, dim = p, mu = mu, Omega = Omega, alpha = alpha)
      fm <- msnFit(x, trace = FALSE, title = NULL, description = "x")
      center <- as.vector(fm@fit$estimated$beta)
      Scatter <- fm@fit$estimated$Omega
      z <- WH.student(x, center = center, cov = Scatter, eta = 0)
      elims[,i] <- sort(z)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band
  }

  n <- nrow(object@data)
  z <- WH.student(object@data, center = c(object@fit$estimated$beta), cov = object@fit$estimated$Omega, eta = 0)

  if (plot.it) {
    band  <- envel(n, c(object@fit$estimated$beta), object@fit$estimated$Omega, object@fit$estimated$alpha, reps, conf)
    ylim <- range(z, band)
    qqnorm(z, ylim = ylim, main = main, font.main = 1, cex.lab = 1.3)
    par(new = TRUE)
    qqnorm(band[,1], axes = FALSE, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = FALSE, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
  }
  invisible(list(transformed = z, envelope = band))
}

envel.st <- function(object, reps = 50, conf = 0.95, plot.it = TRUE, main = "Transformed distances Q-Q plot")
{ ## simulated envelope
  wilson <- function(D2, p, df) {
    F <- D2 / p
    op <- 2 / (9 * p)
    od <- 2 / (9 * df)
    z <- ((1 - od) * F^(1/3) - (1 - op)) / sqrt(op + od * F^(2/3))
    z
  }
  envel2 <- function(n, mu, Omega, alpha, df, reps, conf) {
    conf <- 1 - conf
    p <- length(mu)
    # initialize progress bar
    cat(" Progress:\n")
    pb <- txtProgressBar(min = 0, max = reps, style = 3)
    elims <- matrix(0, nrow = n, ncol = reps)
    for (i in 1:reps) {
      x <- rmvst(n, dim = p, mu = mu, Omega = Omega, alpha = alpha, df = df)
      fm <- mstFit(x, trace = FALSE, title = NULL, description = "x")
      nu <- fm@fit$estimated$nu
      D2 <- Mahalanobis(x, c(fm@fit$estimated$beta), fm@fit$estimated$Omega, inverted = FALSE)
      z <- wilson(D2, p, nu)
      elims[,i] <- sort(z)
      # update progress bar
      setTxtProgressBar(pb, i)
    }
    close(pb)
    band <- matrix(0, nrow = n, ncol = 2)
    for (i in 1:n)
      band[i,] <- quantile(elims[i,], probs = c(conf / 2, 1 - conf / 2))
    band
  }

  n <- nrow(object@data)
  p <- ncol(object@data)
  D2 <- Mahalanobis(object@data, c(object@fit$estimated$beta), object@fit$estimated$Omega, inverted = FALSE)
  z <- wilson(D2, p, object@fit$estimated$nu)

  if (plot.it) {
    band  <- envel2(n, c(object@fit$estimated$beta), object@fit$estimated$Omega, object@fit$estimated$alpha, object@fit$estimated$nu, reps, conf)
    ylim <- range(z, band)
    qqnorm(z, ylim = ylim, main = main, font.main = 1, cex.lab = 1.3)
    par(new = TRUE)
    qqnorm(band[,1], axes = FALSE, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
    par(new = TRUE)
    qqnorm(band[,2], axes = FALSE, main = "", xlab = "", ylab = "", ylim = ylim, type = "l", lwd = 2, col = "red")
  }
  o <- list(transformed = z)
  if (plot.it)
    o$envelope = band
  invisible(o)
}
