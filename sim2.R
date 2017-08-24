### Simulation 2:
### =============
### Generate survival times from a proportional hazards model with time-dependent covariates
### and two strata, and 
### 1) estimate the correctly specified Cox model
### 2) estimate an unstratified Cox model using the stratification variable as a time-independent
###    covariate

library("survival")
library("ggplot2")
source("tvc.R")

set.seed(1234)

N <- 8000

## generate covariates, the strata, and censoring times
stat <- sim_static(1:N) # Time invariant covariates
stat$s <- sample(c(0, 1), N, replace = TRUE, prob = c(0.75, 0.25)) # stratification variable
cens <- sim_cens(1:N) # Censoring times
step <- sim_tvc_step(1:N, cens) # Time-dependent cov: steps
cont <- sim_tvc_cont(1:N, cens) # Time-dependent cov: continuous
tvc <- merge_tvc(step, cont)

## the hazard function of the second stratum is not 'well-behaved'
h1 <- function(t) {
  ifelse(t <= 3, 0.1, 
         ifelse(t <= 10, 0.06, 
                ifelse(t <= 16, 0.15, 0.04)))
}

## simulate event times
evnt <- sapply(1:N, function(i) {
  st <- stat[stat$id == i, c("x1", "x2")]
  tv <- tvc[tvc$id == i, c("x1", "x2", "x3")]
  t <- tvc$time[tvc$id == i]
  s <- stat[stat$id == i, "s"]
  if (s == 0) {
    h <- hazard_gen(baseline = list(type = "gompertz", lambda = 0.05, alpha = 0.04),
                    static = list(X = st, b = c(0.3, 0.3)),
                    step = list(X = tv, b = c(-0.5, -0.3, 0.6), tstart = t))
  } else {
    h <- hazard_gen(baseline = h1,
                    static = list(X = st, b = c(0.3, 0.3)),
                    step = list(X = tv, b = c(-0.5, -0.3, 0.6), tstart = t))
  }
  rsurvt(h)
})
evnt <- event_table(1:N, cens, evnt)

## tmerge the data
st_data <- tmerge(stat, evnt, id = id, tstop = time, event = event(time, event))
tv_data <- tmerge(st_data, tvc, id = id, x.tv1 = tdc(time, x1), 
                  x.tv2 = tdc(time, x2), x.step = tdc(time, x3))

## 1) estimate the stratified model
fit1 <- coxph(Surv(tstart, tstop, event) ~ strata(s) + x1 + x2 + x.tv1 + x.tv2 + x.step, 
              data = tv_data)

## compare estimated parameters to the originals
data.frame(bhat = fit1$coefficients, se = sqrt(diag(fit1$var)),
           b = c(0.3, 0.3, -0.5, -0.3, 0.6)) %>%
  ggplot(aes(x = attr(fit1$coefficients, "names"), y = b)) + 
  geom_errorbar(aes(ymin=bhat-2*se, ymax=bhat+2*se), width = .05) +
  geom_point(aes(y = bhat), color = "black") + geom_point(aes(y = b), color = "red")

## compare estimated baseline survival curves to the original survival functions
newdat <- data.frame(id = c(1, 2), x1 = c(0, 0), x2 = c(0, 0), tstart = c(0, 0), tstop = c(50, 50), 
                     event = c(0, 0), x.tv1 = c(0, 0), x.tv2 = c(0, 0), x.step = c(0, 0), 
                     s = c(0, 1))
plot(survfit(fit1, newdata = newdat), xmax = 40, conf.int = TRUE)
surv0 <- function(t) exp(0.05 / 0.04 * (1 - exp(0.04 * t)))
surv1 <- survival(h1)
curve(surv0, from = 0, to = 40, add = TRUE, col = "red", lwd = 2, lty = 2)
curve(surv1, from = 0, to = 40, add = TRUE, col = "blue", lwd = 2, lty = 2)
