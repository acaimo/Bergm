
.exchange_engine <- function(theta, free, eta_fn, logprior, sy, sim_one,
                             sim_many = NULL, burn.in, main.iters, thin = 1L,
                             gamma, U_noise, adapt = TRUE, target.accept = 0.25,
                             da = NULL, impute = NULL, checkpoint = NULL,
                             blocked = FALSE, verbose = TRUE, out_fn = identity) {
  nch <- ncol(theta)
  p <- nrow(theta)
  pf <- sum(free)
  q <- nrow(sy)
  thin <- max(1L, as.integer(thin))
  n_keep <- main.iters %/% thin
  Theta <- array(NA_real_, c(n_keep, p, nch), dimnames = list(NULL, rownames(theta), NULL))
  if (is.null(sim_many)) sim_many <- function(hs, etas) Map(sim_one, hs, etas)

  eta <- matrix(0, q, nch)
  lp <- numeric(nch)
  for (h in seq_len(nch)) {
    eta[, h] <- eta_fn(theta[, h])
    lp[h] <- logprior(theta[, h])
  }
  if (any(!is.finite(lp)))
    stop("Some starting values have zero prior density; supply 'startVals'.", call. = FALSE)

  log_gamma <- log(gamma)
  acc <- integer(nch)
  n_prop <- 0L
  n_sim <- 0L
  n_sim_burn <- 0L
  n_da_reject <- 0L
  n_fail <- 0L
  use_da <- FALSE
  sur <- NULL
  sur_cur <- numeric(nch)
  da_store <- NULL
  if (!is.null(da)) {
    if (burn.in < 10L) {
      warning("Delayed acceptance needs a burn-in of at least 10 iterations to fit ",
              "its surrogate; it has been switched off.", call. = FALSE)
      da <- NULL
    } else {
      da_start <- floor(burn.in / 2) + 1L
      da_store <- matrix(NA_real_, (burn.in - da_start + 1L) * nch, pf)
      da_row <- 0L
    }
  }

  tot <- burn.in + main.iters
  pb <- .progress(verbose, tot)
  on.exit(pb$close(), add = TRUE)
  groups <- if (blocked) list(seq_len(nch %/% 2), (nch %/% 2 + 1L):nch) else NULL

  finish <- function(h, prop, ep, lpp, s1, s) {
    if (is.null(s)) {
      n_fail <<- n_fail + 1L
      return(FALSE)
    }
    logr <- .exch_ratio(ep, eta[, h], sy[, h], s) + lpp - lp[h] - s1
    if (log(stats::runif(1)) < logr) {
      theta[, h] <<- prop
      eta[, h] <<- ep
      lp[h] <<- lpp
      if (use_da) sur_cur[h] <<- sur(prop)
      TRUE
    } else FALSE
  }

  for (k in seq_len(tot)) {
    in_burn <- k <= burn.in
    g <- exp(log_gamma)
    E <- crossprod(U_noise, matrix(stats::rnorm(pf * nch), pf, nch))
    acc_iter <- 0L

    if (!blocked) {
      for (h in seq_len(nch)) {
        ab <- sample.int(nch - 1L, 2L)
        ab <- ab + (ab >= h)
        prop <- theta[, h]
        prop[free] <- prop[free] + g * (theta[free, ab[1L]] - theta[free, ab[2L]]) + E[, h]
        lpp <- logprior(prop)
        ok <- FALSE
        n_prop <- n_prop + 1L
        if (is.finite(lpp)) {
          s1 <- 0
          go <- TRUE
          if (use_da) {
            s1 <- sur(prop) - sur_cur[h]
            go <- log(stats::runif(1)) < s1
          }
          if (go) {
            ep <- eta_fn(prop)
            if (!anyNA(ep)) {
              s <- sim_one(h, ep)
              if (in_burn) n_sim_burn <- n_sim_burn + 1L else n_sim <- n_sim + 1L
              ok <- finish(h, prop, ep, lpp, s1, s)
            }
          } else n_da_reject <- n_da_reject + 1L
        }
        if (ok) {
          acc_iter <- acc_iter + 1L
          if (!in_burn) acc[h] <- acc[h] + 1L
        }
        if (!is.null(impute)) sy[, h] <- impute(h, theta[, h], ok)
      }
    } else {
      for (gi in 1:2) {
        grp <- groups[[gi]]
        comp <- groups[[3L - gi]]
        props <- vector("list", length(grp))
        lpps <- s1s <- numeric(length(grp))
        need <- logical(length(grp))
        for (j in seq_along(grp)) {
          h <- grp[j]
          ab <- comp[sample.int(length(comp), 2L)]
          prop <- theta[, h]
          prop[free] <- prop[free] + g * (theta[free, ab[1L]] - theta[free, ab[2L]]) + E[, h]
          props[[j]] <- prop
          lpps[j] <- logprior(prop)
          n_prop <- n_prop + 1L
          if (!is.finite(lpps[j])) next
          if (use_da) {
            s1s[j] <- sur(prop) - sur_cur[h]
            if (log(stats::runif(1)) >= s1s[j]) {
              n_da_reject <- n_da_reject + 1L
              next
            }
          }
          need[j] <- TRUE
        }
        eps_all <- vector("list", length(grp))
        for (j in which(need)) {
          eps_all[[j]] <- eta_fn(props[[j]])
          if (anyNA(eps_all[[j]])) need[j] <- FALSE
        }
        idx <- which(need)
        okv <- logical(length(grp))
        if (length(idx)) {
          eps <- eps_all[idx]
          ss <- sim_many(grp[idx], eps)
          if (in_burn) n_sim_burn <- n_sim_burn + length(idx) else n_sim <- n_sim + length(idx)
          for (r in seq_along(idx)) {
            j <- idx[r]
            h <- grp[j]
            if (finish(h, props[[j]], eps[[r]], lpps[j], s1s[j], ss[[r]])) {
              okv[j] <- TRUE
              acc_iter <- acc_iter + 1L
              if (!in_burn) acc[h] <- acc[h] + 1L
            }
          }
        }
        if (!is.null(impute)) for (j in seq_along(grp)) sy[, grp[j]] <- impute(grp[j], theta[, grp[j]], okv[j])
      }
    }

    if (in_burn) {
      if (adapt) {
        log_gamma <- log_gamma + (acc_iter / nch - target.accept) / (k + 1)^0.6
        log_gamma <- min(max(log_gamma, log(1e-3)), log(10))
      }
      if (!is.null(da) && k >= da_start) {
        da_store[da_row + seq_len(nch), ] <- t(theta[free, , drop = FALSE])
        da_row <- da_row + nch
      }
      if (!is.null(da) && k == burn.in) {
        S <- stats::cov(da_store) * da$inflate
        if (all(is.finite(S)) && .is_pd(S)) {
          sur_mean <- colMeans(da_store)
          sur <- .gauss_surrogate(sur_mean, chol(S), free)
          for (h in seq_len(nch)) sur_cur[h] <- sur(theta[, h])
          use_da <- TRUE
          da$mean <- sur_mean
          da$cov <- S
        } else {
          warning("Delayed-acceptance surrogate could not be fitted (degenerate burn-in); ",
                  "continuing without it.", call. = FALSE)
        }
        da_store <- NULL
      }
    } else {
      kk <- k - burn.in
      if (kk %% thin == 0L) Theta[kk %/% thin, , ] <- out_fn(theta)
    }

    if (!is.null(checkpoint) && k %in% checkpoint$at) {
      partialBergmEstimate <- c(list(Theta = Theta, failedAfter = k), checkpoint$extra())
      save(partialBergmEstimate, file = checkpoint$file)
    }
    pb$tick(k)
  }

  list(Theta = Theta, acc = acc, gamma = exp(log_gamma), n_prop = n_prop,
       n_sim = n_sim, n_sim_burn = n_sim_burn, n_da_reject = n_da_reject, n_fail = n_fail,
       da = if (use_da) da else NULL, theta = out_fn(theta), sy = sy)
}

.gauss_surrogate <- function(mean, U, free) {
  force(mean); force(U); force(free)
  function(th) {
    r <- backsolve(U, th[free] - mean, transpose = TRUE)
    -0.5 * sum(r * r)
  }
}

.bergm_output <- function(res, setup, offset, prior, main.iters, thin, cut.reject,
                          formula, runtime, extra = list()) {
  if (isTRUE(res$n_fail > 0L))
    warning(res$n_fail, " proposal(s) were rejected because the auxiliary network exceeded ",
            "an ergm limit (edge count or term cut-off); the model may be near-degenerate ",
            "in parts of the parameter space.", call. = FALSE)
  Theta <- res$Theta
  nch <- dim(Theta)[3]
  names_ <- setup$coef.names
  dimnames(Theta)[[2]] <- names_
  chains <- coda::mcmc.list(lapply(seq_len(nch), function(h)
    coda::mcmc(matrix(Theta[, , h], ncol = length(names_), dimnames = list(NULL, names_)),
               thin = thin)))
  pooled <- do.call(rbind, lapply(chains, as.matrix))
  colnames(pooled) <- names_
  if (cut.reject) pooled <- unique(pooled)
  FF <- coda::mcmc(pooled)

  free <- !offset
  ess <- rep(NA_real_, length(names_))
  names(ess) <- names_
  ess[free] <- round(coda::effectiveSize(chains[, free, drop = FALSE]), 0)

  rhat <- rep(NA_real_, length(names_))
  names(rhat) <- names_
  if (nch >= 2L && nrow(Theta) >= 4L) {
    rhat[free] <- tryCatch(
      coda::gelman.diag(chains[, free, drop = FALSE], autoburnin = FALSE,
                        multivariate = FALSE)$psrf[, 1],
      error = function(e) NA_real_)
  }
  fixed <- offset
  names(fixed) <- names_
  AR <- round(sum(res$acc) / (main.iters * nch), 2)
  AR.chains <- round(res$acc / main.iters, 3)

  out <- c(list(Time = runtime,
                formula = formula,
                specs = names_,
                dim = length(names_),
                Theta = FF,
                AR = AR,
                ess = ess,
                fixed = fixed,
                chains = chains,
                Rhat = rhat,
                AR.chains = AR.chains,
                nchains = nch,
                gamma = res$gamma,
                prior = list(mean = prior$mean, sigma = prior$sigma, type = prior$type,
                             decay = prior$decay),
                constraints = setup$constraints,
                n.sims = c(burn.in = res$n_sim_burn, main = res$n_sim),
                n.failed = res$n_fail,
                delayed.acceptance = if (!is.null(res$da))
                  list(stage1.rejections = res$n_da_reject, mean = res$da$mean,
                       cov = res$da$cov) else NULL),
           extra)
  class(out) <- "bergm"
  out
}
