
.pl_setup <- function(formula, constraints = NULL) {
  m <- suppressMessages(ergmMPLE(formula, constraints = constraints %||% ~., output = "matrix"))
  X <- m$predictor
  storage.mode(X) <- "double"
  list(X = X, y = as.numeric(m$response), w = as.numeric(m$weights),
       names = colnames(X), p = ncol(X))
}

.pl_loglik <- function(pl, theta) cpp_logpl(pl$X, pl$y, pl$w, as.numeric(theta))

.pl_grad <- function(pl, theta) {
  eta <- drop(pl$X %*% theta)
  drop(crossprod(pl$X, pl$w * (pl$y - stats::plogis(eta))))
}

.pl_hess <- function(pl, theta) {
  pr <- stats::plogis(drop(pl$X %*% theta))
  H <- -crossprod(pl$X, pl$X * (pl$w * pr * (1 - pr)))
  dimnames(H) <- list(pl$names, pl$names)
  H
}

.pl_mode <- function(pl, init = rep(0, pl$p), prior = NULL, maxit = 200L, tol = 1e-10) {
  P <- if (is.null(prior)) matrix(0, pl$p, pl$p) else prior$prec
  mu <- if (is.null(prior)) rep(0, pl$p) else prior$mean
  obj <- function(th) .pl_loglik(pl, th) - 0.5 * sum((th - mu) * drop(P %*% (th - mu)))
  th <- as.numeric(init)
  f <- obj(th)
  for (it in seq_len(maxit)) {
    g <- .pl_grad(pl, th) - drop(P %*% (th - mu))
    H <- .pl_hess(pl, th) - P
    step <- tryCatch(solve(-H, g), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) return(NULL)
    a <- 1
    repeat {
      th_new <- th + a * step
      f_new <- obj(th_new)
      if (is.finite(f_new) && f_new >= f - 1e-12) break
      a <- a / 2
      if (a < 1e-10) break
    }
    if (max(abs(th_new)) > 1e4) return(NULL)
    conv <- abs(f_new - f) < tol * (abs(f) + tol) && max(abs(a * step)) < 1e-8
    th <- th_new
    f <- f_new
    if (conv) break
  }
  names(th) <- pl$names
  th
}

.gauss_prior <- function(prior.mean, prior.sigma, p, names) {
  pr <- .make_prior(prior.mean, prior.sigma, p, rep(TRUE, p), names)
  list(mean = unname(pr$mean), sigma = unname(pr$sigma), prec = pr$prec, U = pr$U,
       logdens = function(Theta) .dmvn_log(Theta, unname(pr$mean), pr$U),
       full = pr)
}

.rwm_pl <- function(pl, prior, S.prop, init, burn.in, main.iters, thin = 1L,
                    temperature = 1, adj = NULL) {
  p <- pl$p
  W <- adj$W %||% diag(p)
  tm <- adj$Theta_MLE %||% rep(0, p)
  tp <- adj$Theta_PL %||% rep(0, p)
  L <- t(.chol_pd(S.prop))
  out <- cpp_rwm_adjpl(pl$X, pl$y, pl$w, unname(as.matrix(W)), as.numeric(tm), as.numeric(tp),
                       as.numeric(temperature), prior$mean, prior$prec, as.numeric(init),
                       L, as.integer(burn.in), as.integer(main.iters), as.integer(thin))
  colnames(out$draws) <- pl$names
  out
}

.adjpl_batch <- function(pl, Theta, adj = NULL, grad = FALSE) {
  p <- pl$p
  if (is.null(dim(Theta))) Theta <- matrix(Theta, nrow = 1)
  W <- adj$W %||% diag(p)
  cpp_adjpl_batch(pl$X, pl$y, pl$w, unname(as.matrix(Theta)), unname(as.matrix(W)),
                  as.numeric(adj$Theta_MLE %||% rep(0, p)),
                  as.numeric(adj$Theta_PL %||% rep(0, p)), grad)
}

.map_newton <- function(st, sy, init, prior, aux.iters, n.aux.draws, aux.thin,
                        iters = 20L, verbose = FALSE) {
  p <- length(init)
  sample_at <- function(th, state, burnin) {
    z <- .simulate_state(state, th, burnin = burnin, samplesize = n.aux.draws,
                         interval = aux.thin)
    if (is.null(z)) return(NULL)
    S <- z$s
    g <- sy - colMeans(S) - drop(prior$prec %*% (th - prior$mean))
    H <- .make_pd(stats::cov(S) + prior$prec)
    md <- sqrt(max(0, sum(g * solve(H, g))))
    list(th = th, g = g, H = H, md = md, state = update(state, state = z$state))
  }
  cur <- sample_at(as.numeric(init), st, aux.iters)
  if (is.null(cur)) stop(.sim_limit_msg, call. = FALSE)
  rewarm <- max(aux.thin, aux.iters %/% 4L)
  radius <- 1
  traj <- matrix(NA_real_, iters, p)
  accepted <- logical(iters)
  half <- max(1L, iters %/% 2L)
  for (i in seq_len(iters)) {
    step <- solve(cur$H, cur$g)
    md_step <- sqrt(sum(step * drop(cur$H %*% step)))
    gain <- if (i <= half) 1 else 1 / (i - half + 1)
    a <- gain * min(1, radius / max(md_step, 1e-12))
    new <- sample_at(cur$th + a * step, cur$state, rewarm)
    ok <- !is.null(new) && is.finite(new$md) && new$md <= cur$md + 0.5
    if (ok) {
      cur <- new
      radius <- min(2, radius * 1.5)
    } else {
      radius <- radius / 2
    }
    accepted[i] <- ok
    traj[i, ] <- cur$th
  }
  keep <- (half + 1L):iters
  list(theta = colMeans(traj[keep, , drop = FALSE]), trajectory = traj,
       state = cur$state, accepted = accepted, discrepancy = cur$md)
}

.map_rm <- function(st, sy, init, prior, aux.iters, n.aux.draws, aux.thin,
                    rm.iters, rm.a, rm.alpha) {
  p <- length(init)
  theta <- matrix(0, rm.iters, p)
  theta[1, ] <- init
  for (i in 2:rm.iters) {
    z <- .simulate_state(st, theta[i - 1, ], burnin = aux.iters,
                         samplesize = n.aux.draws, interval = aux.thin)
    if (is.null(z)) stop(.sim_limit_msg, call. = FALSE)
    estgrad <- -colMeans(sweep(z$s, 2, sy)) - drop(prior$prec %*% (theta[i - 1, ] - prior$mean))
    theta[i, ] <- theta[i - 1, ] + (rm.a / i) * (rm.alpha + estgrad)
  }
  list(theta = theta[rm.iters, ], trajectory = theta, state = st)
}

.ergm_estimate <- function(formula, estimate, seed, constraints = NULL, ...) {
  fit <- suppressMessages(ergm(formula, estimate = estimate, verbose = FALSE,
                               constraints = constraints %||% ~.,
                               control = control.ergm(seed = seed), ...))
  stats::coef(fit)
}
