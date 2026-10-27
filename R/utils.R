
`%||%` <- function(a, b) if (is.null(a)) b else a

.as_cov <- function(S, p, what = "covariance") {
  if (is.null(S)) return(NULL)
  if (is.matrix(S)) {
    if (!all(dim(S) == c(p, p)))
      stop(sprintf("'%s' must be a %d x %d matrix.", what, p, p), call. = FALSE)
    return(unname(S))
  }
  if (length(S) == 1L) return(diag(as.numeric(S), p, p))
  if (length(S) == p) return(diag(as.numeric(S), p, p))
  stop(sprintf("'%s' must be a scalar, a vector of length %d or a %d x %d matrix.",
               what, p, p, p), call. = FALSE)
}

.make_pd <- function(S, tol = 1e-8) {
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  floor <- max(tol, tol * max(abs(e$values)))
  if (all(e$values > floor)) return(S)
  v <- pmax(e$values, floor)
  out <- e$vectors %*% (v * t(e$vectors))
  dimnames(out) <- dimnames(S)
  (out + t(out)) / 2
}

.is_pd <- function(S, tol = 1e-10) {
  if (!all(is.finite(S))) return(FALSE)
  ev <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
  all(ev > tol * max(1, max(abs(ev))))
}

.chol_pd <- function(S) {
  out <- tryCatch(chol(S), error = function(e) NULL)
  if (is.null(out)) out <- chol(.make_pd(S))
  out
}

.rmvn <- function(n, mean, U) {
  p <- length(mean)
  Z <- matrix(stats::rnorm(n * p), n, p)
  sweep(Z %*% U, 2, mean, "+")
}

.dmvn_log <- function(x, mean, U) {
  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  p <- ncol(x)
  R <- backsolve(U, t(x) - mean, transpose = TRUE)
  -0.5 * colSums(R * R) - sum(log(diag(U))) - 0.5 * p * log(2 * pi)
}

.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.logmeanexp <- function(x) .logsumexp(x) - log(length(x))

.make_prior <- function(prior.mean, prior.sigma, p, free, coef.names,
                        log.prior = NULL, dec = rep(FALSE, p), prior.decay = NULL) {
  nd <- sum(dec)
  expand_mean <- function(m) {
    if (nd == 0L || length(m) != p - nd) return(m)
    out <- rep(0, p)
    out[!dec] <- m
    out
  }
  expand_sigma <- function(S) {
    if (nd == 0L || is.null(S)) return(S)
    k <- p - nd
    if (is.matrix(S) && all(dim(S) == k)) {
      out <- diag(1, p)
      out[!dec, !dec] <- S
      return(out)
    }
    if (!is.matrix(S) && length(S) == k) {
      out <- rep(1, p)
      out[!dec] <- S
      return(out)
    }
    S
  }
  prior.mean <- if (is.null(prior.mean)) rep(0, p) else expand_mean(as.numeric(prior.mean))
  if (length(prior.mean) == 1L && p > 1L) prior.mean <- rep(prior.mean, p)
  if (length(prior.mean) != p)
    stop(sprintf("'prior.mean' must have length %d (the number of model parameters).", p),
         call. = FALSE)
  prior.sigma <- .as_cov(expand_sigma(prior.sigma %||% 100), p, "prior.sigma")
  g <- free & !dec
  m <- prior.mean[g]
  S <- prior.sigma[g, g, drop = FALSE]
  if (length(m) && !.is_pd(S)) stop("'prior.sigma' must be positive definite.", call. = FALSE)
  U <- if (length(m)) chol(S) else matrix(0, 0, 0)
  P <- if (length(m)) chol2inv(U) else matrix(0, 0, 0)
  const <- -sum(log(diag(U))) - 0.5 * length(m) * log(2 * pi)
  names(prior.mean) <- coef.names
  dimnames(prior.sigma) <- list(coef.names, coef.names)

  decf <- free & dec
  pd <- .decay_prior(prior.decay, sum(decf), coef.names[decf])

  if (!is.null(log.prior)) {
    if (!is.function(log.prior)) stop("'log.prior' must be a function.", call. = FALSE)
    logdens <- function(theta) {
      jac <- 0
      if (any(decf)) {
        jac <- sum(theta[decf])
        theta[decf] <- exp(theta[decf])
      }
      names(theta) <- coef.names
      v <- as.numeric(log.prior(theta))
      if (length(v) != 1L || is.na(v)) -Inf else v + jac
    }
    type <- "custom"
  } else {
    logdens <- function(theta) {
      out <- 0
      if (length(m)) {
        r <- backsolve(U, theta[g] - m, transpose = TRUE)
        out <- const - 0.5 * sum(r * r)
      }
      if (any(decf))
        out <- out + sum(stats::dnorm(theta[decf], pd$meanlog, pd$sdlog, log = TRUE))
      out
    }
    type <- "normal"
  }
  if (any(dec)) {
    prior.mean[dec] <- NA_real_
    prior.sigma[dec, ] <- prior.sigma[, dec] <- NA_real_
  }
  list(mean = prior.mean, sigma = prior.sigma, mean.free = m, sigma.free = S,
       U = U, prec = P, logdens = logdens, type = type,
       decay = if (any(decf)) pd)
}

.decay_prior <- function(prior.decay, k, names) {
  pd <- prior.decay %||% c(meanlog = log(0.5), sdlog = 0.75)
  pd <- as.list(pd)
  if (is.null(pd$meanlog) || is.null(pd$sdlog))
    stop("'prior.decay' must give 'meanlog' and 'sdlog'.", call. = FALSE)
  ml <- as.numeric(pd$meanlog)
  sl <- as.numeric(pd$sdlog)
  if (k > 0L) {
    if (!length(ml) %in% c(1L, k) || !length(sl) %in% c(1L, k))
      stop(sprintf("'prior.decay' entries must have length 1 or %d.", k), call. = FALSE)
    if (any(sl <= 0)) stop("'prior.decay$sdlog' must be positive.", call. = FALSE)
    ml <- rep_len(ml, k)
    sl <- rep_len(sl, k)
    names(ml) <- names(sl) <- names
  }
  list(meanlog = ml, sdlog = sl)
}

.decay_index <- function(setup, offset = setup$offset) {
  dec <- rep(FALSE, setup$p)
  for (cv in setup$etamap$curved) {
    idx <- cv$from[-1L]
    idx <- idx[grepl("\\.(decay|lambda)$", setup$coef.names[idx])]
    dec[idx] <- TRUE
  }
  dec & !offset
}

.param_map <- function(setup, dec) {
  eta <- .eta_fun(setup)
  to_nat <- function(x) {
    if (is.matrix(x)) x[dec, ] <- exp(x[dec, , drop = FALSE]) else x[dec] <- exp(x[dec])
    x
  }
  to_work <- function(x) {
    if (is.matrix(x)) x[dec, ] <- log(x[dec, , drop = FALSE]) else x[dec] <- log(x[dec])
    x
  }
  list(dec = dec, to_nat = to_nat, to_work = to_work,
       eta = if (any(dec)) function(x) eta(to_nat(x)) else eta)
}

.ergm_setup <- function(formula, constraints = NULL, control.args = list(),
                        network = NULL) {
  if (!is.null(network)) formula <- .formula_with_network(formula, network)
  y <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  p <- nparam(model, canonical = FALSE)
  if (is.null(constraints)) constraints <- ~.
  control.args <- .term_opts(control.args, network.size(y))
  ctrl <- do.call(control.simulate.formula, c(list(MCMC.burnin = 1L, MCMC.interval = 1L),
                                              control.args))
  sim <- suppressMessages(
    simulate(formula, coef = rep(0, p), nsim = 1, constraints = constraints,
             control = ctrl, return.args = "ergm_state"))
  state <- sim$object
  model <- state$model
  coef.names <- param_names(model, canonical = FALSE)
  list(formula = formula, network = y, state = state, model = model,
       etamap = model$etamap, p = p, q = nparam(model, canonical = TRUE),
       coef.names = coef.names, offset = model$etamap$offsettheta,
       curved = is.curved(model), sy = as.numeric(state$stats),
       stat.names = param_names(model, canonical = TRUE),
       constraints = constraints, directed = is.directed(y),
       n = network.size(y))
}

.eta_fun <- function(setup) {
  if (!setup$curved) return(function(theta) theta)
  etamap <- setup$etamap
  function(theta) ergm.eta(theta, etamap)
}

.SLAVE_CONTROL <- list(MCMC.maxedges = Inf, MCMC.save_networks = FALSE)

.simulate_state <- function(st, eta, burnin, samplesize = 1L, interval = 1L) {
  z <- tryCatch(
    ergm_MCMC_slave(st, eta = eta, control = .SLAVE_CONTROL, verbose = FALSE,
                    burnin = as.integer(burnin), samplesize = as.integer(samplesize),
                    interval = as.integer(interval)),
    error = function(e) {
      if (grepl("cut-off|cutoff", conditionMessage(e))) NULL else stop(e)
    })
  if (is.null(z)) return(NULL)
  if (!identical(as.integer(z$status), 0L)) {
    if (identical(as.integer(z$status), 1L)) return(NULL)
    stop("Network simulation failed (Metropolis-Hastings proposal failure).",
         call. = FALSE)
  }
  z
}

.exch_ratio <- function(eta1, eta0, sy, s) {
  d <- eta1 - eta0
  ok <- is.finite(d) & d != 0
  if (!any(ok)) return(0)
  sum(d[ok] * (sy[ok] - s[ok]))
}

.check_offsets <- function(setup, offset.coef) {
  off <- setup$offset
  if (!is.null(offset.coef)) {
    if (any(is.na(offset.coef))) stop("NaN, NA are not allowed in offset.coef.", call. = FALSE)
    if (length(offset.coef) != sum(off))
      stop("Invalid offset parameter vector offset.coef: wrong number of parameters: expected ",
           sum(off), " got ", length(offset.coef), ".", call. = FALSE)
  } else if (any(off)) {
    stop("The model contains offset terms whose parameter values have not been specified: ",
         paste(setup$coef.names[off], collapse = ", "), ".", call. = FALSE)
  }
  invisible(TRUE)
}

.mple <- function(formula, setup, offset.coef = NULL, constraints = NULL) {
  fit <- tryCatch(
    suppressMessages(suppressWarnings(
      ergm(formula, estimate = "MPLE", verbose = FALSE, offset.coef = offset.coef,
           constraints = constraints %||% ~.))),
    error = function(e) NULL)
  p <- setup$p
  if (is.null(fit)) {
    return(list(coef = stats::setNames(rep(0, p), setup$coef.names), vcov = NULL,
                ok = FALSE))
  }
  cf <- stats::coef(fit)
  V <- tryCatch(as.matrix(stats::vcov(fit)), error = function(e) NULL)
  list(coef = cf, vcov = V, ok = TRUE)
}

.fix_infinite <- function(mple, offset, offset.coef, verbose = TRUE) {
  if (!any(is.infinite(mple))) return(list(offset = offset, offset.coef = offset.coef))
  vals <- rep(NA_real_, length(offset))
  if (!is.null(offset.coef)) vals[offset] <- offset.coef
  free_idx <- which(!offset)
  ninf <- which(mple[free_idx] == -Inf)
  pinf <- which(mple[free_idx] == Inf)
  if (length(ninf) && verbose)
    message("Observed statistic(s) for ", paste(names(mple)[free_idx][ninf], collapse = ", "),
            " are at their smallest attainable values; their coefficients are fixed at -Inf.")
  if (length(pinf) && verbose)
    message("Observed statistic(s) for ", paste(names(mple)[free_idx][pinf], collapse = ", "),
            " are at their largest attainable values; their coefficients are fixed at Inf.")
  vals[free_idx[ninf]] <- -Inf
  vals[free_idx[pinf]] <- Inf
  offset[free_idx[c(ninf, pinf)]] <- TRUE
  list(offset = offset, offset.coef = as.numeric(vals[offset]))
}

.bergm_worker_env <- new.env(parent = emptyenv())

.worker_set <- function(name, value) {
  assign(name, value, envir = .bergm_worker_env)
  NULL
}

.worker_sim <- function(etas, burnin, samplesize = 1L, interval = 1L, key = "state") {
  st <- get(key, envir = .bergm_worker_env)
  lapply(etas, function(eta) {
    z <- .simulate_state(st, eta, burnin, samplesize, interval)
    if (is.null(z)) NULL else z$s
  })
}

.worker_sim_states <- function(jobs, burnin, samplesize = 1L, interval = 1L) {
  lapply(jobs, function(j) {
    z <- .simulate_state(j$state, j$eta, burnin, samplesize, interval)
    if (is.null(z)) NULL else z$s
  })
}

.start_cluster <- function(cores, seed = NULL) {
  cl <- parallel::makePSOCKcluster(cores)
  this_lib <- dirname(getNamespaceInfo("Bergm", "path"))
  parallel::clusterCall(cl, ".libPaths", unique(c(this_lib, .libPaths())))
  parallel::clusterEvalQ(cl, suppressPackageStartupMessages({
    requireNamespace("ergm"); requireNamespace("Bergm")
  }))
  parallel::clusterSetRNGStream(cl, seed)
  cl
}

.chunk <- function(x, k) {
  k <- max(1L, min(k, length(x)))
  split(x, cut(seq_along(x), k, labels = FALSE))
}

.par_lapply <- function(X, FUN, ..., cl = NULL) {
  if (is.null(cl) || length(X) <= 1L) return(lapply(X, FUN, ...))
  parallel::clusterApplyLB(cl, X, FUN, ...)
}

.with_seed <- function(seed) {
  if (!is.null(seed) && !is.na(seed)) set.seed(seed)
  invisible(NULL)
}

.progress <- function(verbose, total) {
  if (!isTRUE(verbose) || !interactive()) return(list(tick = function(k) NULL, close = function() NULL))
  pb <- utils::txtProgressBar(min = 0, max = total, style = 3)
  list(tick = function(k) utils::setTxtProgressBar(pb, k), close = function() close(pb))
}

.term_opts <- function(control.args, n) {
  to <- control.args$term.options %||% list()
  if (is.null(to$gw.cutoff)) to$gw.cutoff <- max(30L, min(as.integer(n) - 2L, 1000L))
  control.args$term.options <- to
  control.args
}

.formula_rhs <- function(formula) paste(deparse(formula[[3]], width.cutoff = 500L), collapse = " ")

.formula_with_network <- function(formula, network) {
  f <- formula
  f[[2]] <- as.name(".bergm_nw")
  environment(f) <- list2env(list(.bergm_nw = network), parent = environment(formula))
  f
}

.not_supported_curved <- function(setup, fun) {
  if (setup$curved)
    stop(fun, "() does not support curved ERGM terms; use fixed decay parameters ",
         "(e.g. gwesp(0.5, fixed = TRUE)) or bergm().", call. = FALSE)
  if (any(setup$offset))
    stop(fun, "() does not support offset terms; use bergm().", call. = FALSE)
}
