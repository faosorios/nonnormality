# loading required packages
library(fastmatrix)
library(plot3D)

# reading R sources
source("mardia.R")

# Figure 11.a
par(pty = "s", mai = c(0.5,0.5,0.5,0.5))
persp3D(z = z, facets = FALSE, theta = 115, phi = 20, xlab = "gamma", ylab = "nu", zlab = "skewness", bty = "b2", lwd = 2, colkey = FALSE)

# Figure 11.b
par(pty = "s", mai = c(0.5,0.5,0.5,0.5))
persp3D(z = z, facets = FALSE, theta = 210, phi = 20, xlab = "gamma", ylab = "nu", zlab = "skewness", bty = "b2", lwd = 2, colkey = FALSE)
