# loading required packages
library(MVT)
library(fMultivar)

# loading dataset and reading R sources
data(WindSpeed)
source("../code/boots.R")
source("../code/envelope.R")
source("../code/stats.R")
source("../code/student.influence.R")

# save dimensions
n <- nrow(WindSpeed)
p <- ncol(WindSpeed)

## Table 2: Descriptive statistics
stats(WindSpeed)
#$summary.1D
#       mean       sd   skewness kurtosis
#vs 16.98009 13.61276 -0.8492518 4.375535
#gh 12.74052 13.33339 -0.6919762 2.617095
#kw 14.03223 17.23853 -0.4096161 3.008603
#
#$Scatter
#         vs       gh       kw
#vs 185.3073 126.9589 148.1797
#gh 126.9589 177.7792 110.6200
#kw 148.1797 110.6200 297.1668
#
#$mardia
#[1]  3.516002 23.849807

# Gaussian fit
fn <- studentFit(WindSpeed, family = Student(eta = 0))

## Table 3
fn
#Call:
#studentFit(x = WindSpeed, family = Student(eta = 0))
#Converged in 1 iterations
#
#Center:
#      vs      gh      kw 
#16.9801 12.7405 14.0322 
#
#Scatter matrix estimate:
#   vs       gh       kw      
#vs 184.6408                  
#gh 126.5022 177.1397         
#kw 147.6466 110.2221 296.0979
#
#Number of Observations: 278 

# log-likelihood
fn$logLik
#[1] -3254.534

# Multivariate t fit
ft <- studentFit(WindSpeed, family = Student(eta = .25))

## Table 4 (estimate of eta is reported in 'Call' output)
ft
#Call:
#studentFit(x = WindSpeed, family = Student(eta = 0.2365))
#Converged in 28 iterations
#
#Center:
#      vs      gh      kw 
#18.9587 14.8311 16.7253 
#
#Scatter matrix estimate:
#   vs       gh       kw      
#vs 217.2405                  
#gh 175.1545 246.0688         
#kw 200.7992 152.0299 353.5857
#
#Number of Observations: 278 

# log-likelihood
ft$logLik
#[1] -3211.186

## Figure 1.a
obs <- c(16,42,47,233)
D2 <- fn$distances / p
par(pty = "s")
plot(D2, ylim = c(0,12), ylab = "Modified Mahalanobis distances", cex.lab = 1.3)
cutoff <- qgamma(p = 0.975, shape = p / 2, rate = p / 2)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text((1:n)[obs], D2[obs], label = as.character(obs), pos = 3)

## Figure 1.b.
par(pty = "s")
envelope.student(fn, reps = 5000)

# Figure 2.a
eta <- ft$eta
F <- (ft$distances / p) / (1 - 2 * eta)
cutoff <- qf(p = 0.975, df1 = p, df2 = 1 / eta)
par(pty = "s")
plot(F, ylim = c(0,25), ylab = "Modified Mahalanobis distances", cex.lab = 1.3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")
text((1:n)[obs], F[obs], label = as.character(obs), pos = 3)

## Figure 2.b.
par(pty = "s")
envelope.student(ft, reps = 5000)

# bootstrap confidence intervals
o <- boots(WindSpeed)

## Table 5 (results may vary slightly because I forgot to assign the seed, my bad)
o[c("eta", "neg", "ci.eta", "ci.neg")]
#$eta
#[1] 0.2365211
#
#$neg
#[1] 0.2787819
#
#$ci.eta
#                   L         U
#normal     0.1661919 0.3068503
#percentile 0.1639751 0.3043134
#pivotal    0.1687287 0.3090670
#
#$ci.neg
#                    L         U
#normal     0.06533974 0.4922241
#percentile 0.11906954 0.5349825
#pivotal    0.02258141 0.4384943

## Figure 3.a.
par(pty = "s")
hist(o$stat[,1], main = NULL, breaks = 20, col = "gray", xlab = "eta", freq = FALSE, cex.lab = 1.3)
xgrid <- seq(from = 0, to = 0.5, len = 300)
ygrid <- dnorm(xgrid, mean = o$eta, sd = sd(o$stat[,1]))
lines(xgrid, ygrid, col = "red", lwd = 2)

## Figure 3.b.
par(pty = "s")
hist(o$stat[,2], main = NULL, breaks = 20, col = "gray", xlab = "Negentropy", freq = FALSE, cex.lab = 1.3)
xgrid <- seq(from = 0, to = 1, len = 300)
ygrid <- dnorm(xgrid, mean = o$neg, sd = sd(o$stat[,2]))
lines(xgrid, ygrid, col = "red", lwd = 2)

# Multivariate t fit: observations 16,42,47,233 removed
ss <- (1:n)[-obs]
f0 <- studentFit(WindSpeed, family = Student(eta = .25), subset = ss)

## Table 6
f0
#Call:
#studentFit(x = WindSpeed, family = Student(eta = 0.1892), subset = ss)
#Converged in 60 iterations
#
#Center:
#      vs      gh      kw 
#18.8326 14.6204 16.6128 
#
#Scatter matrix estimate:
#   vs       gh       kw      
#vs 188.6560                  
#gh 154.9214 217.4122         
#kw 176.5383 132.3037 310.0284
#
#Number of Observations: 274 

# log-likelihood
f0$logLik
#[1] -3135.325

# Local influence: t-perturbation 
z <- influence.t(as.matrix(WindSpeed))

## Figure 4.a
par(pty = "s")
plot(z)

## Figure 4.b
obs <- c(16,42,233)
cutoff <- mean(z$total) + 2 * sd(z$total)
par(pty = "s")
plot(z$total, ylim = c(0,1), ylab = "diag(B)", cex.lab = 1.3)
text((1:n)[obs], z$total[obs], as.character(obs), pos = 3)
abline(h = cutoff, lwd = 2, lty = 2, col = "red")

# skew-normal fit
fsn <- msnFit(WindSpeed, description = "skew normal")

## Table 7
fsn
#
#Title:
# Skew Normal Parameter Estimation 
#
#Call:
# msnFit(x = WindSpeed, description = "skew normal")
#
#Model:
# Skew Normal Distribution
#
#Estimated Parameter(s):
#$beta
#           vs       gh       kw
#[1,] 25.79198 28.20022 23.11244
#
#$Omega
#         vs       gh       kw
#vs 262.2901 262.7314 227.6605
#gh 262.7314 416.1421 250.5995
#kw 227.6605 250.5995 378.5482
#
#$alpha
#        vs         gh         kw 
# 1.0728846 -4.9740819 -0.2818751 

# log-likelihood
fsn@fit$logL
#[1] -3229.176

# skew-t fit
fst <- mstFit(WindSpeed, description = "skew t")

## Table 8
fst
#
#Title:
# Student-t Parameter Estimation 
#
#Call:
# mstFit(x = WindSpeed, description = "skew t")
#
#Model:
# Skew Student-t Distribution
#
#Estimated Parameter(s):
#$beta
#          vs      gh       kw
#[1,] 27.4857 27.8276 24.32632
#
#$Omega
#         vs       gh       kw
#vs 175.5073 188.0342 161.1409
#gh 188.0342 276.3197 166.3621
#kw 161.1409 166.3621 236.5383
#
#$alpha
#        vs         gh         kw 
# 0.6328237 -4.5479609 -0.1238838 
#
#$nu
#[1] 4.046509

# log-likelihood
fst@fit$logL
#[1] -3180.724

## loading required package
library(sn)

# bootstrap confidence intervals, skew-normal model
o <- boots.st(WindSpeed)

## Table 9: columns 1 and 2
o[c("neg","CI")]
#$neg
#[1] 0.2459538
#
#$CI
#                   L         U
#normal     0.1882251 0.3036826
#percentile 0.1763347 0.3027681
#pivotal    0.1891396 0.3155729

# bootstrap confidence intervals, skew-t model (lengthy operation)
o <- boots.st(WindSpeed)

## Table 9: columns 3 to 6
o[c("nu","neg","CI.nu","CI.neg")]
#$nu
#[1] 4.046509
#
#$neg
#[1] 0.4410002
#
#$CI.nu
#                  L        U
#normal     2.559610 5.533408
#percentile 3.239191 5.946032
#pivotal    2.146986 4.853827
#
#$CI.neg
#                   L         U
#normal     0.3310159 0.5509845
#percentile 0.3265290 0.5419484
#pivotal    0.3400520 0.5554714

## Figure 5.a
par(pty = "s")
envel.sn(fsn, reps = 5000)

## Figure 5.b (this is a lengthy operation)
par(pty = "s")
envel.st(fst, reps = 5000)
