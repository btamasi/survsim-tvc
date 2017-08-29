### Simulation 4:
### =============
## TODO: write description

library("survival")
library("ggplot2")
source("tvc.R")

set.seed(98765)

N <- 5000
p <- 0.2

## generate covariates and censoring times
stat <- sim_static(1:N) # Time invariant covariates
cens <- sim_cens(1:N) # Censoring times
tvc <- sim_tvc_cont(1:N, cens) # Time-dependent cov: continuous

## the hazard function is not 'well-behaved'
h0 <- function(t) {
  ifelse(t <= 4, 0.1, 
         ifelse(t <= 12, 0.06,
                ifelse(t <= 20, 0.02, 
                       ifelse(t <= 30, 0.15, 0.08))))
}

## simulate event times
evnt <- sapply(1:N, function(i) {
  st <- stat[stat$id == i, c("x1", "x2")]
  tv <- tvc[tvc$id == i, c("x1", "x2")]
  t <- tvc$time[tvc$id == i]
  h <- hazard_gen(baseline = h0,
                    static = list(X = st, b = c(0.4, -0.3)),
                    step = list(X = tv, b = c(-0.3, 0.6), tstart = t))
  rsurvt(h)
})
evnt <- event_table(1:N, cens, evnt)

## choose observations that countinue to live
id_s <- sample(evnt$id[evnt$event == 0], round(sum(evnt$event == 0) * p))
id_s <- id_s[order(id_s)]
stat_s <- stat[stat$id %in% id_s, ] # Time invariant covariates
cens_s <- sim_cens(id_s) # Censoring times
tvc_s <- sim_tvc_cont(id_s, cens_s) # Time-dependent cov: continuous

## simulate event times
evnt_s <- sapply(id_s, function(i) {
  st <- stat_s[stat_s$id == i, c("x1", "x2")]
  tv <- tvc_s[tvc_s$id == i, c("x1", "x2")]
  t <- tvc_s$time[tvc_s$id == i]
  h <- hazard_gen(baseline = list(type = "exp", lambda = 0.02),
                  static = list(X = st, b = c(0.6, -0.15)),
                  step = list(X = tv, b = c(-0.5, 0.3), tstart = t))
  rsurvt(h)
})
evnt_s <- event_table(id_s, cens_s, evnt_s)

evnt_s$time <- evnt$time[evnt$id %in% id_s] + evnt_s$time
tvc_s$time <- do.call(c, lapply(id_s, function(i) {
  tvc_s$time[tvc_s$id == i] + evnt$time[evnt$id == i]
}))

evnt[evnt$id %in% id_s, ] <- evnt_s
tvc$s <- rep(0, nrow(tvc))
tvc_s$s <- rep(1, nrow(tvc_s))
tvc <- rbind(tvc, tvc_s)
tvc <- tvc[order(tvc$id, tvc$time), ]

## tmerge the data
st_data <- tmerge(stat, evnt, id = id, tstop = time, event = event(time, event))
tv_data <- tmerge(st_data, tvc, id = id, x.tv1 = tdc(time, x1), 
                  x.tv2 = tdc(time, x2), s = tdc(time, s))

## estimate the misspecified model
fit <- coxph(Surv(tstart, tstop, event) ~ x1 + x2 + x.tv1 + x.tv2 + s, 
              data = tv_data)
summary(fit)
