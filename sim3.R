### Simulation 3:
### =============
### Generate survival times from a proportional hazards model with time-dependent covariates
### and a set of redundant covariates. Estimate the model.

library("survival")
library("ggplot2")
source("tvc.R")

set.seed(13)

N <- 8000

## generate covariates and censoring times
stat <- sim_static(1:N) # Time invariant covariates
stat_r <- sim_static(1:N, varnames = c("xr1", "xr2")) # Time invariant redundant covariates
cens <- sim_cens(1:N) # Censoring times
step <- sim_tvc_step(1:N, cens) # Time-dependent cov: steps
step_r <- sim_tvc_step(1:N, cens, varname = "xr3") # Redundant time-dependent cov: steps
cont <- sim_tvc_cont(1:N, cens) # Time-dependent cov: continuous
cont_r <- sim_tvc_cont(1:N, cens, varnames = c("xr4", "xr5")) # Redundant time-dependent cov: continuous
tvc <- merge_tvc(step, cont)
tvc_r <- merge_tvc(step_r, cont_r)

## simulate event times
evnt <- sapply(1:N, function(i) {
  st <- stat[stat$id == i, c("x1", "x2")]
  tv <- tvc[tvc$id == i, c("x1", "x2", "x3")]
  t <- tvc$time[tvc$id == i]
  h <- hazard_gen(baseline = list(type = "weibull", lambda = 0.06, nu = 1.5),
                  static = list(X = st, b = c(0.3, 0.3)),
                  step = list(X = tv, b = c(-0.5, -0.3, 0.6), tstart = t))
  rsurvt(h)
})
evnt <- event_table(1:N, cens, evnt)

## tmerge the data
st_data <- tmerge(merge(stat, stat_r, all = TRUE), 
                  evnt, id = id, tstop = time, event = event(time, event))
tv_data <- tmerge(st_data, tvc, id = id, x.tv1 = tdc(time, x1), 
                  x.tv2 = tdc(time, x2), x.step = tdc(time, x3))
tv_data_r <- tmerge(tv_data, tvc_r, id = id, xr.tv1 = tdc(time, xr4), 
                    xr.tv2 = tdc(time, xr5), xr.step = tdc(time, xr3))

## fit the model
fit <- coxph(Surv(tstart, tstop, event) ~ x1 + x2 + x.tv1 + x.tv2 + x.step + 
               xr1 + xr2 + xr.tv1 + xr.tv2 + xr.step, data = tv_data_r)

## compare estimated parameters to the originals
data.frame(bhat = fit$coefficients, se = sqrt(diag(fit$var)),
           b = c(0.3, 0.3, -0.5, -0.3, 0.6, 0, 0, 0, 0, 0)) %>%
  ggplot(aes(x = attr(fit$coefficients, "names"), y = b)) + 
  geom_errorbar(aes(ymin=bhat-2*se, ymax=bhat+2*se), width = .05) +
  geom_point(aes(y = bhat), color = "black") + geom_point(aes(y = b), color = "red")

## compare estimated baseline survival curve to the original survival function
newdat <- data.frame(id = 1, x1 = 0, x2 = 0, tstart = 0, tstop = 50, event = 0, 
                     x.tv1 = 0, x.tv2 = 0, x.step = 0, xr1 = 0, xr2 = 0, 
                     xr.tv1 = 0, xr.tv2 = 0, xr.step = 0)
plot(survfit(fit, newdata = newdat))
curve(exp(-0.06*x^1.5), from = 0, to = 50, add = TRUE, col = "red", lty = 2, lwd = 2) # Weibull survival curve
