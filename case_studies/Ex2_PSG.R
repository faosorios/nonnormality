# loading required packages
library(MVT)
library(fMultivar)

# loading dataset and reading R sources
data(PSG)
source("../code/boots.R")
source("../code/envelope.R")
source("../code/stats.R")
source("../code/student.influence.R")

# Gaussian fit
fn <- studentFit(~ manual + automated, data = PSG, family = Student(eta = 0))

# extract 'x' matrix and save dimensions
x <- fn$x
n <- nrow(x)
p <- ncol(x)

## Table 10: Descriptive statistics
stats(x)
#$summary.1D
#              mean        sd skewness kurtosis
#manual    2.553891 0.8781231 1.527482 2.770493
#automated 2.308982 1.1190175 5.889866 3.120458
#
#$Scatter
#             manual automated
#manual    0.7711002 0.7027743
#automated 0.7027743 1.2522002
#
#$mardia
#[1]  1.005507 18.291620

## Table 11
fn
#Call:
#studentFit(x = ~manual + automated, data = PSG, family = Student(eta = 0))
#Converged in 1 iterations
#
#Center:
#    manual automated 
#   2.5539    2.3090 
#
#Scatter matrix estimate:
#          manual    automated
#manual    0.7616965          
#automated 0.6942039 1.2369294
#
#Number of Observations: 82 

# log-likelihood
fn$logLik
#[1] -200.8901

# Multivariate t fit
ft <- studentFit(~ manual + automated, data = PSG, family = Student(eta = 0.25))

## Table 12 (estimate of eta is reported in 'Call' output)
ft
#Call:
#studentFit(x = ~manual + automated, data = PSG, family = Student(eta = 0.4528))
#Converged in 132 iterations
#
#Center:
#    manual automated 
#   2.6149    2.5299 
#
#Scatter matrix estimate:
#          manual   automated
#manual    4.825547          
#automated 4.749887 5.217942 
#
#Number of Observations: 82 

# log-likelihood
ft$logLik
#[1] -165.7818

# bootstrap confidence intervals
o <- boots(x)

## Table 13 (results may vary slightly because I forgot to assign the seed, my bad)
o[c("eta", "neg", "ci.eta", "ci.neg")]
#$eta
#[1] 0.4527691
#
#$neg
#[1] 1.454022
#
#$ci.eta
#                   L         U
#normal     0.3631895 0.5423487
#percentile 0.3057744 0.4896735
#pivotal    0.4158647 0.5997638
#
#$ci.neg
#                     L        U
#normal     0.193868883 2.714174
#percentile 0.334048523 2.900547
#pivotal    0.007495976 2.573994

## Figure 6.a
obs <- c(1,30,79)
D2 <- fn$distances / p
par(pty = "s")
plot(D2, ylim = c(0,12), ylab = "Modified Mahalanobis distances", cex.lab = 1.3)
cutoff <- qgamma(p = 0.975, shape = p / 2, rate = p / 2)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text((1:n)[obs], D2[obs], label = as.character(obs), pos = 3)

## Figure 6.b.
par(pty = "s")
envelope.student(fn, reps = 5000)

# Figure 7.a
obs <- c(1,30,35,79)
eta <- ft$eta
F <- (ft$distances / p) / (1 - 2 * eta)
cutoff <- qf(p = 0.975, df1 = p, df2 = 1 / eta)
par(pty = "s")
plot(F, ylim = c(0,140), ylab = "Modified Mahalanobis distances", cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text((1:n)[obs], F[obs], label = as.character(obs), pos = 3)

## Figure 7.b.
par(pty = "s")
envelope.student(ft, reps = 5000)

## Figure 8.a.
par(pty = "s")
hist(o$stat[,1], main = NULL, breaks = 20, col = "gray", xlab = "eta", freq = FALSE, cex.lab = 1.3)
xgrid <- seq(from = 0, to = 0.5, len = 300)
ygrid <- dnorm(xgrid, mean = o$eta, sd = sd(o$stat[,1]))
lines(xgrid, ygrid, col = "red", lwd = 2)

## Figure 8.b.
par(pty = "s")
hist(o$stat[,2], main = NULL, breaks = 20, col = "gray", xlab = "Negentropy", freq = FALSE, cex.lab = 1.3)
xgrid <- seq(from = 0, to = 4, len = 300)
ygrid <- dnorm(xgrid, mean = o$neg, sd = sd(o$stat[,2]))
lines(xgrid, ygrid, col = "red", lwd = 2)

# Local influence: t-perturbation 
z <- influence.t(x)

## Figure 9.a
par(pty = "s")
plot(z)

## Figure 9.b
obs <- c(30,79)
cutoff <- mean(z$total) + 2 * sd(z$total)
par(pty = "s")
plot(z$total, ylim = c(0,1), ylab = "diag(B)", cex.lab = 1.3)
text((1:n)[obs], z$total[obs], as.character(obs), pos = 3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")

# skew-normal fit
fsn <- msnFit(x, description = "skew normal")

## Table 14
fsn
#
#Title:
# Skew Normal Parameter Estimation 
#
#Call:
# msnFit(x = x, description = "skew normal")
#
#Model:
# Skew Normal Distribution
#
#Estimated Parameter(s):
#$beta
#       manual automated
#[1,] 2.246923  1.696943
#
#$Omega
#             manual automated
#manual    0.8559259 0.8820803
#automated 0.8820803 1.6115209
#
#$alpha
#    manual  automated 
#-0.1100308  0.8457029 

# log-likelihood
fsn@fit$logL
#[1] -200.7823

# skew-t fit
fst <- mstFit(x, description = "skew t")

## Table 15
fst
#
#Title:
# Student-t Parameter Estimation 
#
#Call:
# mstFit(x = x, description = "skew t")
#
#Model:
# Skew Student-t Distribution
#
#Estimated Parameter(s):
#$beta
#       manual automated
#[1,] 2.822276  2.929567
#
#$Omega
#             manual automated
#manual    0.4647986 0.4853370
#automated 0.4853370 0.5843961
#
#$alpha
#    manual  automated 
#  8.436844 -10.423684 
#
#$nu
#[1] 1.749065

# log-likelihood
fst@fit$logL
#[1] -151.4756

## loading required package
library(sn)

# bootstrap confidence intervals, skew-normal model
o <- boots.sn(x)

## Table 16: columns 1 and 2
o[c("neg","CI")]
#$neg
#[1] 0.05020276
#
#$CI
#                    L          U
#normal     -0.2657429 0.36614847
#percentile  0.0248734 0.50137202
#pivotal    -0.4009665 0.07553213

# bootstrap confidence intervals, skew-t model (lengthy operation)
o <- boots.st(x)

## Table 16: columns 3 to 6
o[c("nu","neg","CI.nu","CI.neg")]
#$nu
#[1] 1.749065
#
#$neg
#[1] 0.5940079
#
#$CI.nu
#                   L        U
#percentile  1.153714 4.526732
#pivotal    -1.028601 2.344417
#
#$CI.neg
#                   L         U
#percentile 0.3624674 0.8396852
#pivotal    0.3483306 0.8255484

## Figure 10.a
par(pty = "s")
envel.sn(fsn, reps = 5000)

## Figure 10.b (this is a lengthy operation)
par(pty = "s")
envel.st(fst, reps = 5000)
