sim_static <- function(ids, varnames = NULL) {
  N <- length(ids)
  x1 <- runif(N, -0.5, 0.5)
  x2 <- sample(c(0, 1), N, replace = TRUE, prob = c(0.3, 0.7))
  val <- data.frame(id = ids, x1, x2)
  if (!is.null(varnames))
    names(val) <- c("id", varnames)
  return(val)
}


sim_cens <- function(ids, ec = 25) {
  ## Generate exponentially distributed censoring times for each id
  N <- length(ids)
  rexp(N, 1/ec)
}


sim_tvc_step <- function(ids, ct, er = 15, varname = NULL) {
  N <- length(ids)
  rt <- rexp(N, 1/er)
  idr <- which(ct > rt)
  id <- rep(ids, (ct > rt) + 1)
  x3 <- c(rep(0, N), rep(1, length(idr)))
  x3 <- x3[order(c(1:N, idr + 0.5))]
  t <- c(rep(0, N), rt[idr])
  t <- t[order(c(1:N, idr + 0.5))]
  val <- data.frame(id, time = t, x3)
  if (!is.null(varname))
    names(val) <- c("id", "time", varname)
  return(val)
}


sim_tvc_cont <- function(ids, ct, varnames = NULL) {
  N <- length(ids)
  id <- rep(ids, floor(ct)+1)
  t <- lapply(ct, function(x) seq(0, floor(x)))
  t <- do.call(c, t)
  x1 <- rnorm(length(t))
  x2 <- sample(c(0, 1), length(t), replace = TRUE, prob = c(0.8, 0.2))
  val <- data.frame(id, time = t, x1, x2)
  if (!is.null(varnames))
    names(val) <- c("id", "time", varnames)
  return(val)
}


merge_tvc <- function(...) {
  merge_all <- function(x, y) merge(x, y, all = TRUE)
  val <- Reduce(merge_all, list(...))
  val <- sapply(val, function(x) {
    good_id <- !is.na(x)
    good_val <- x[good_id]
    good_val[cumsum(good_id)]
  })
  data.frame(val)
}


event_table <- function(id, cens, survt) {
  data.frame(id, time = pmin(cens, survt), event = as.numeric(survt < cens))
}


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


cum_hazard <- function(hazard, subdiv = 2000) {
  ## TODO: an option to tweak with tolerance
  function(t) {
    integrate(hazard, 0, t, subdivisions = subdiv)$value
  }
}


survival <- function(hazard) {
  function(t) {
    ch <- cum_hazard(hazard)
    val <- sapply(t, ch)
    exp(-val)
  }
}


rsurvt <- function(hazard, upper = 1000) {
  u <- runif(1)
  ch <- cum_hazard(hazard)
  fun <- function(t) { -log(u) - ch(t) }
  uniroot(fun, interval = c(0, upper))$root
}


survsim <- function(x, beta, h0 = list(type = "exp", lambda = 0.1), 
                    cens = NULL, discrete = FALSE) {
  x <- as.matrix(x)
  N <- nrow(x)
  stopifnot(ncol(x) == length(beta))
  ## Draw survival times
  rs <- exp(x %*% beta)
  U <- runif(N)
  Tsim <- switch(h0$type,
                 exp = -log(U) / (h0$lambda * rs),
                 weibull = (-log(U) / (h0$lambda * rs)) ^ (1 / h0$nu),
                 gompertz = 1 / h0$alpha * log(1 - (h0$alpha * log(U)) / (h0$lambda * rs))
  )
  ## Censoring
  if (!is.null(cens)) {
    evnt <- as.numeric(Tsim <= cens)
    Tsim <- pmin(Tsim, cens)
  } else {
    evnt <- rep(1, N)
  }
  ## Discrete times
  if (discrete) {
    Tsim <- ceiling(Tsim)
  }
  data.frame(time = Tsim, 
             event = evnt,
             x)
}