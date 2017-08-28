### Simulation 2/B:
### ===============
### Generate survival times from a proportional hazards model with time-invariant covariates
### and two strata, and estimate an unstratified Cox model using the stratification variable
### as a covariate. Replicate the simulation many times to see what problems this misspecification
### can cause. The time invariant version is used for speeding up the simulations.

library("survival")
library("ggplot2")
source("tvc.R")


## 1) Estimate the model and take a look at baseline survivals
set.seed(111222)

## Simulate the data
N <- 5000
X <- data.frame(x1 = rnorm(N, 3, 2), x2 = sample(c(0, 1), N, replace = TRUE, prob = c(0.4, 0.6)),
                s = sample(c(0, 1), N, replace = TRUE, prob = c(0.25, 0.75)))
b <- c(0.3, -0.5)
C <- rexp(N, 1/30) # Exponentially distributed censoring times
T1 <- survsim(X[X$s == 0, c("x1", "x2")], b, h0 = list(type = "exp", lambda = 0.03), 
              cens = C[X$s == 0])
T1$s <- rep(0, nrow(T1))
T2 <- survsim(X[X$s == 1, c("x1", "x2")], b, 
              h0 = list(type = "gompertz", lambda = 0.05, alpha = 0.04), 
              cens = C[X$s == 1])
T2$s <- rep(1, nrow(T2))
simdata <- rbind(T1, T2)

## Fit the stratified model
fit <- coxph(Surv(time, event) ~ strata(s) + x1 + x2, data = simdata)
summary(fit)

## Plot the survival curves
newdat <- data.frame(x1 = c(0, 0), x2 = c(0, 0), s = c(0, 1))
plot(survfit(fit, newdata = newdat), xmax = 40, conf.int = TRUE)
surv1 <- function(t) exp(-0.03 * t)
surv2 <- function(t) exp(0.05 / 0.04 * (1 - exp(0.04 * t)))
curve(surv1, from = 0, to = 40, add = TRUE, col = "red", lwd = 2, lty = 2)
curve(surv2, from = 0, to = 40, add = TRUE, col = "blue", lwd = 2, lty = 2)

## 2) Simulate the misspecified model
set.seed(112233)

N <- 10000
reps <- 5000
b <- c(0.3, -0.5)

par <- replicate(reps, {
  X <- data.frame(x1 = rnorm(N, 3, 2), x2 = sample(c(0, 1), N, replace = TRUE, prob = c(0.4, 0.6)),
                  s = sample(c(0, 1), N, replace = TRUE, prob = c(0.25, 0.75)))
  C <- rexp(N, 1/30) # Exponentially distributed censoring times
  T1 <- survsim(X[X$s == 0, c("x1", "x2")], b, h0 = list(type = "exp", lambda = 0.01), 
                cens = C[X$s == 0])
  T1$s <- rep(0, nrow(T1))
  T2 <- survsim(X[X$s == 1, c("x1", "x2")], b, 
                h0 = list(type = "gompertz", lambda = 0.05, alpha = 0.04), 
                cens = C[X$s == 1])
  T2$s <- rep(1, nrow(T2))
  simdata <- rbind(T1, T2)
  fit <- coxph(Surv(time, event) ~ s + x1 + x2, data = simdata)
  coef(fit)[c("x1", "x2")]
})

par(mfrow = c(1, 2))
hist(par["x1", ], breaks = 40, col = grey(0.6, alpha = 0.75), lty = 0,
     panel.first = grid(col = "black"), main = "Histogram of b1", xlab = "b1")
rug(par["x1", ], col = grey(0.6, alpha = 0.2))
abline(v = b[1], col = "red", lwd = 2)

hist(par["x2", ], breaks = 40, col = grey(0.6, alpha = 0.75), lty = 0,
     panel.first = grid(col = "black"), main = "Histogram of b2", xlab = "b2")
rug(par["x2", ], col = grey(0.6, alpha = 0.2))
abline(v = b[2], col = "red", lwd = 2)

## 3) Is the problem more pronounced when the stratifcation variable is correlated with
## one of the covariates? -- No clear evidence
set.seed(998877)

N <- 10000
reps <- 5000
b <- c(0.3, -0.5)
rho <- -0.9 # correlation coefficient

par <- replicate(reps, {
  x1 <- rnorm(N)
  s <- rho * x1 + sqrt(1 - rho^2) * rnorm(N)
  x1 <- 3 + 2 * x1
  s <- as.numeric(pnorm(s) > 0.25)
  X <- data.frame(x1 = x1, x2 = sample(c(0, 1), N, replace = TRUE, prob = c(0.4, 0.6)),
                  s = s)
  C <- rexp(N, 1/30) # Exponentially distributed censoring times
  T1 <- survsim(X[X$s == 0, c("x1", "x2")], b, h0 = list(type = "exp", lambda = 0.01), 
                cens = C[X$s == 0])
  T1$s <- rep(0, nrow(T1))
  T2 <- survsim(X[X$s == 1, c("x1", "x2")], b, 
                h0 = list(type = "gompertz", lambda = 0.05, alpha = 0.04), 
                cens = C[X$s == 1])
  T2$s <- rep(1, nrow(T2))
  simdata <- rbind(T1, T2)
  fit <- coxph(Surv(time, event) ~ s + x1 + x2, data = simdata)
  coef(fit)[c("x1", "x2")]
})

par(mfrow = c(1, 2))
hist(par["x1", ], breaks = 40, col = grey(0.6, alpha = 0.75), lty = 0,
     panel.first = grid(col = "black"), main = "Histogram of b1", xlab = "b1")
rug(par["x1", ], col = grey(0.6, alpha = 0.2))
abline(v = b[1], col = "red", lwd = 2)

hist(par["x2", ], breaks = 40, col = grey(0.6, alpha = 0.75), lty = 0,
     panel.first = grid(col = "black"), main = "Histogram of b2", xlab = "b2")
rug(par["x2", ], col = grey(0.6, alpha = 0.2))
abline(v = b[2], col = "red", lwd = 2)
