sim_covs <- function(N, ec = 50, er = 15) {
  id <- NULL
  tstart <- NULL
  tstop <- NULL
  x1 <- NULL
  x2 <- NULL
  x3 <- NULL
  for (i in 1:N) {
    ## Censoring time
    ct <- rexp(1, 1/ec)
    ## Regime change time
    rt <- rexp(1, 1/er)
    if (ct > rt) {
      ## ids
      id <- c(id, rep(i, 2))
      ## time
      tstart <- c(tstart, c(0, rt))
      tstop <- c(tstop, c(rt, ct))
      ## static covairates
      x1 <- c(x1, rep(runif(1, -0.5, 0.5), 2))
      x2 <- c(x2, rep(sample(c(0, 1), 1, prob = c(0.7, 0.3)), 2))
      ## time dependent covariate
      x3 <- c(x3, c(0, 1))
    } else {
      ## ids
      id <- c(id, i)
      ## time
      tstart <- c(tstart, 0)
      tstop <- c(tstop, ct)
      ## static covairates
      x1 <- c(x1, runif(1, -0.5, 0.5))
      x2 <- c(x2, sample(c(0, 1), 1, prob = c(0.7, 0.3)))
      ## time dependent covariate
      x3 <- c(x3, 0)
    }
  }
  data.frame(id, tstart, tstop, x1, x2, x3)
}

sim_static <- function(ids) {
  N <- length(ids)
  x1 <- runif(N, -0.5, 0.5)
  x2 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.3, 0.7))
  data.frame(id = ids, x1, x2)
}

sim_cens <- function(ids, ec = 25) {
  ## Generate exponentially distributed censoring times for each id
  N <- length(ids)
  rexp(N, 1/ec)
}

sim_tvc_step <- function(ids, ct, er = 15) {
  N <- length(ids)
  rt <- rexp(N, 1/er)
  idr <- which(ct > rt)
  id <- rep(ids, (ct > rt) + 1)
  x3 <- c(rep(0, N), rep(1, length(idr)))
  x3 <- x3[order(c(1:N, idr + 0.5))]
  t <- c(rep(0, N), rt[idr])
  t <- t[order(c(1:N, idr + 0.5))]
  data.frame(id, time = t, x3)
}


sim_tvc_cont <- function(ids, ct) {
  N <- length(ids)
  id <- rep(ids, floor(ct)+1)
  t <- lapply(ct, function(x) seq(0, floor(x)))
  t <- do.call(c, t)
  x1 <- rnorm(length(t))
  x2 <- sample(c(0, 1), length(t), replace = TRUE, prob = c(0.8, 0.2))
  data.frame(id, time = t, x1, x2)
}

merge_tvc <- function(x, y) {
  val <- merge(x, y, all = TRUE)
  val <- sapply(val, function(x) {
    good_id <- !is.na(x)
    good_val <- x[good_id]
    good_val[cumsum(good_id)]
  })
  data.frame(val)
}

sim_event <- function(static, tvc) {
  ids <- static$id
  et <- numeric(length(ids))
  for (i in seq_along(ids)) {
    id <- ids[i]
    haz <- hazard_gen(baseline = list(type = "exp", lambda = 0.1), 
                      static = list(X = static[static$id == id, c("x1", "x2")], b = c(0.3, 0.3)),
                      step = list(X = tvc[tvc$id == id, c("x1", "x2", "x3")], b = c(-0.5, -0.3, 0.6), 
                                  tstart = tvc$time[tvc$id == id]))
    et[i] <- rsurvt(haz)
  }
  return(et)
}

event_table <- function(id, cens, survt) {
  data.frame(id, time = pmin(cens, survt), event = as.numeric(survt < cens))
}

tvc_table <- function(tvc, event) {
  val <- mapply(function(x, y) { x[x$time <= y, ] }, 
                split(tvc, tvc$id), event$time, 
                SIMPLIFY = FALSE, USE.NAMES = FALSE)
  do.call(rbind, val)
}

## TODO: rewrite the hazard function generator functions: 
## - baseline hazard: exp, wei, gom, arbitrary function of time
## - time invariant parts
## - step function parts (as a function of time)
## - time varying parts (additionally to the relative hazard e.g. trends
##   and interaction terms for time varying effects of covariates)

hazard_gen <- function(baseline = list(type = "exp", lambda = 0.1),
                       static = NULL,
                       step = NULL,
                       tfun = NULL) {
  ### Return an arbitrary (proportional) hazard function
  if (is.list(baseline)) {
    h0 <- switch(baseline$type,
                 exp = function(t) rep(baseline$lambda, length(t)),
                 weibull = function(t) baseline$lambda * baseline$nu * t ^ (baseline$nu - 1),
                 gompertz = function(t) baseline$lambda * exp(baseline$alpha * t)
                 )
  } else if (is.function(baseline)) {
    h0 <- baseline
  } else if (is.null(baseline)){
    h0 <- function(t) rep(1, length(t))
  } else {
    stop("baseline is either a list or a function")
  }
  if (!is.null(static)) {
    rs <- exp(sum(static$X * static$b))
  } else {
    rs <- 1
  }
  if (!is.null(step)) {
    Xs <- as.matrix(step$X)
    b <- cbind(step$b)
    hstep <- function(t) {
      i <- findInterval(t, step$tstart)
      as.numeric(exp(Xs[i, ] %*% b))
    }
  } else {
    hstep <- function(t) rep(1, length(t))
  }
  if (!is.null(tfun)) {
    htfun <- tfun
  } else {
    htfun <- function(t) rep(1, length(t))
  }
  haz <- function(t) {
    h0(t) * rs * hstep(t) * htfun(t)
  }
  return(haz)
}

haz_weibull_tvc <- function(lambda, nu, X, b, start, stop) {
  X <- as.matrix(X)
  haz <- function(t) {
    if (t >= max(stop)) {
      sc <- exp(sum(X[nrow(X), ] * b))
    } else {
      for (i in seq_along(start)) {
        if (t >= start[i] & t < stop[i]) {
          sc <- exp(sum(X[i, ] * b))
          break
        } 
      }
    }
    lambda * sc * nu * t ^ (nu - 1)
  }
  function(t) sapply(t, haz)
}

haz_exp_tvc <- function(lambda, X, b, start, stop) {
  X <- as.matrix(X)
  haz <- function(t) {
    if (t >= max(stop)) {
      sc <- exp(sum(X[nrow(X), ] * b))
    } else {
      for (i in seq_along(start)) {
        if (t >= start[i] & t < stop[i]) {
          sc <- exp(sum(X[i, ] * b))
          break
        } 
      }
    }
    lambda * sc
  }
  function(t) sapply(t, haz)
}

haz_exp_tvc2 <- function(lambda, x, b, step_time, step_size = 0, trend = 0) {
  function(t) {
    lambda * exp(sum(x*b) + ifelse(t >= step_time, step_size, 0) + trend*t)
  }
}

haz_exp <- function(lambda) {
  function(t) {
    rep(lambda, length(t))
  }
}

Haz <- function(h) {
  function(t) {
    integrate(h, 0, t)$value
  }
}

Haz2 <- function(h) {
  function(t) {
    sapply(t, function(x) integrate(h, 0, x)$value)
  }
}

# Surv <- function(H) {
#   function(t) {
#     exp(-H(t))
#   }
# }

rsurvt2 <- function(H) {
  u <- runif(1)
  fu <- function(t) { -log(u) - H(t) }
  uniroot(fu, interval = c(0, 1000))$root
}

cum_hazard <- function(hazard, subdiv = 2000) {
  ## TODO: an option to tweak with tolerance
  function(t) {
    integrate(hazard, 0, t, subdivisions = subdiv)$value
  }
}

survival <- function(hazard) {
  ## TODO: vectorize
  function(t) {
    ch <- cum_hazard(hazard)
    exp(-ch(t))
  }
}

rsurvt <- function(hazard, upper = 1000) {
  u <- runif(1)
  ch <- cum_hazard(hazard)
  fun <- function(t) { -log(u) - ch(t) }
  uniroot(fun, interval = c(0, upper))$root
}

sim_gen <- function(N, ...) {
  fnc <- deparse(match.call())
  fnc <- sub("^sim_gen\\(N = [0-9]+, ", "", fnc)
  fnc <- sub("\\)$", "", fnc)
  fnc <- sub("\\), ", "\\),,,", fnc)
  fnc <- strsplit(fnc, ",,,")
  fnc <- lapply(fnc, function(x) strsplit(x, " = "))
  n <- length(fnc[[1]])
  vnames <- vector("character", n)
  simdat <- data.frame(matrix(0, nrow = N, ncol = n))
  #browser()
  for (i in 1:n) {
    vnames[i] <- fnc[[1]][[i]][1]
    c <- fnc[[1]][[i]][2]
    vn <- sub("^.+\\([ ]*", "", c)
    vn <- sub("[ ]*,.+)|)$", "", vn)
    ll <- list(N)
    names(ll) <- vn
    simdat[, i] <- eval(parse(text = c), envir = ll)
  }
  names(simdat) <- vnames
  
  # vname <- function(c) {
  #   c <- deparse(f)
  #   vn <- sub("^.+\\([ ]*", "", c)
  #   vn <- sub("[ ]*,.+)|)$", "", vn)
  #   vn
  # }
  
  # fns <- list(...)
  # fns <- lapply(fns, substitute)
  # vns <- sapply(fns, vname)
  # ll <- list(N)
  # names(ll) <- vn
  # eval(f, envir = ll)
  simdat
}
