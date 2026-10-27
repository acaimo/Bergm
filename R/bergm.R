#' Parameter estimation for Bayesian ERGMs
#'
#' Fits Bayesian exponential random graph models with the approximate exchange
#' algorithm and a population of interacting chains updated by adaptive
#' direction sampling (ADS) (Caimo and Friel, 2011).
#'
#' Compared with Bergm 5, the sampler draws every auxiliary network through a
#' single call to \pkg{ergm}'s compiled sampler on a pre-built network state,
#' supports curved ERGM terms, models with a single parameter, custom priors,
#' burn-in adaptation of the ADS move factor, delayed acceptance and parallel
#' chains. The posterior sample is returned both pooled (\code{Theta}, as in
#' previous versions) and per chain (\code{chains}), which enables standard
#' convergence diagnostics such as the potential scale reduction factor.
#'
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}. Curved terms (e.g.
#' \code{gwesp(fixed = FALSE)}) are supported.
#'
#' @param prior.mean vector; mean vector of the multivariate Normal prior.
#' By default set to a vector of 0's.
#'
#' @param prior.sigma square matrix, vector or scalar; variance/covariance
#' matrix for the multivariate Normal prior. A scalar or vector is interpreted
#' as the diagonal. By default a diagonal matrix with every diagonal entry
#' equal to 100.
#'
#' @param burn.in count; number of burn-in iterations for every chain of the
#' population.
#'
#' @param main.iters count; number of iterations for every chain of the
#' population.
#'
#' @param aux.iters count; number of auxiliary iterations used for network
#' simulation.
#'
#' @param nchains count; number of chains of the population MCMC. By default
#' twice the number of free model parameters (at least 4), rounded up to a
#' multiple of \code{2 * cores} when \code{cores > 1}.
#'
#' @param gamma scalar; parallel adaptive direction sampling move factor. When
#' \code{adapt = TRUE} this is the starting value, tuned during burn-in.
#'
#' @param V.proposal scalar, vector or matrix; covariance of the small Normal
#' perturbation added to every ADS move. A scalar is the diagonal entry
#' (default 0.0025, as in previous versions).
#'
#' @param startVals vector; optional starting values for the parameter
#' estimation. Alternatively a previous \code{bergm} fit of the same model, in
#' which case every chain continues from where that run ended (and its adapted
#' \code{gamma} is reused unless \code{gamma} is given), so a run can be
#' extended in stages.
#'
#' @param offset.coef vector; a vector of coefficients for the offset terms.
#'
#' @param constraints formula; constraints on the support of the network
#' distribution, as in \code{\link[ergm]{ergm}}.
#'
#' @param thin count; thinning interval applied to every chain.
#'
#' @param cut.reject logical; if \code{TRUE}, repeated (rejected) draws are
#' removed from the pooled sample \code{Theta}. This changes the posterior
#' summaries and is only kept for backward compatibility; default \code{FALSE}.
#'
#' @param saveEveryX count; if not \code{NULL}, the partial posterior sample is
#' saved to \code{saveEveryXName} every \code{saveEveryX} iterations.
#'
#' @param saveEveryXName character; file name of the partial estimation object.
#'
#' @param seed integer; optional seed for reproducibility (also seeds the
#' parallel workers when \code{cores > 1}).
#'
#' @param adapt logical; tune \code{gamma} during burn-in towards
#' \code{target.accept}. Adaptation stops at the end of burn-in.
#'
#' @param target.accept scalar; target acceptance rate for \code{adapt}.
#'
#' @param delayed.acceptance logical; if \code{TRUE} (default), proposals are first
#' screened with a Gaussian surrogate of the posterior fitted to the second half
#' of the burn-in population, and an auxiliary network is simulated only for
#' proposals that pass. The second stage corrects for the surrogate so the
#' target distribution is unchanged (Christen and Fox, 2005). It needs a
#' burn-in of at least 10 iterations; with a shorter burn-in it is switched off
#' (silently when left at its default). \code{FALSE} gives the Bergm 5 kernel.
#'
#' @param da.inflate scalar; variance inflation factor for the
#' delayed-acceptance surrogate (values > 1 make the screen more permissive).
#'
#' @param cores count; number of parallel workers. With \code{cores > 1} the
#' population is updated in two halves (each half conditionally on the other)
#' and the auxiliary simulations of a half run in parallel. Worthwhile when a
#' single auxiliary simulation takes more than a few milliseconds.
#'
#' @param log.prior function; optional user-defined log prior density taking
#' the full (named) parameter vector, with decay parameters on their natural
#' (positive) scale, and returning a scalar. Overrides the Normal and
#' log-normal priors; return \code{-Inf} outside the support.
#'
#' @param prior.decay named vector or list with elements \code{meanlog} and
#' \code{sdlog}; log-normal prior for the decay parameters of curved terms
#' (e.g. \code{gwesp.decay}, \code{gwdegree.decay}, \code{altkstar.lambda}),
#' recycled over all of them or given one value per decay parameter. Default:
#' median 0.5, \code{sdlog} 0.75. Decay parameters are sampled on the log
#' scale, so they are always positive; \code{prior.mean} and
#' \code{prior.sigma} apply to the other parameters and may be given either
#' for all parameters (decay entries are ignored) or for the non-decay ones only.
#'
#' @param verbose logical; print progress information.
#'
#' @param ... additional arguments passed to
#' \code{\link[ergm]{control.simulate.formula}} to configure the network
#' proposal (e.g. \code{MCMC.prop}, \code{MCMC.prop.weights}).
#'
#' @return An object of class \code{bergm} with elements
#' \item{Theta}{pooled posterior sample (\code{mcmc}), chains stacked;}
#' \item{chains}{the per-chain posterior samples (\code{mcmc.list});}
#' \item{AR, AR.chains}{overall and per-chain acceptance rates;}
#' \item{ess, Rhat}{effective sample sizes and potential scale reduction factors;}
#' \item{fixed}{logical vector flagging fixed (offset) parameters;}
#' \item{gamma}{the (adapted) ADS move factor;}
#' \item{n.sims}{number of auxiliary network simulations performed;}
#' plus the model specification and run settings.
#'
#' @references
#' Caimo, A. and Friel, N. (2011), "Bayesian Inference for Exponential Random Graph Models,"
#' Social Networks, 33(1), 41-55. \url{https://arxiv.org/abs/1007.5192}
#'
#' Caimo, A. and Friel, N. (2014), "Bergm: Bayesian Exponential Random Graphs in R,"
#' Journal of Statistical Software, 61(2), 1-25. \doi{10.18637/jss.v061.i02}
#'
#' Christen, J. A. and Fox, C. (2005), "Markov chain Monte Carlo using an
#' approximation," Journal of Computational and Graphical Statistics, 14(4), 795-810.
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data(florentine)
#'
#' # Posterior parameter estimation:
#' p.flo <- bergm(flomarriage ~ edges + kstar(2),
#'                burn.in    = 50,
#'                aux.iters  = 500,
#'                main.iters = 3000,
#'                gamma      = 1.2,
#'                seed       = 1)
#'
#' # Posterior summaries:
#' summary(p.flo)
#'
#' # Curved model, delayed acceptance and 4 parallel workers:
#' p.laz <- bergm(lazega ~ edges + nodematch("Office") + gwesp(0.5, fixed = FALSE),
#'                burn.in = 200, main.iters = 2000, aux.iters = 3000,
#'                delayed.acceptance = TRUE, cores = 4, seed = 1)
#' }
#' @export
bergm <- function(formula,
                  prior.mean = NULL,
                  prior.sigma = NULL,
                  burn.in = 100,
                  main.iters = 1000,
                  aux.iters = 1000,
                  nchains = NULL,
                  gamma = 0.5,
                  V.proposal = 0.0025,
                  startVals = NULL,
                  offset.coef = NULL,
                  constraints = NULL,
                  thin = 1,
                  cut.reject = FALSE,
                  saveEveryX = NULL,
                  saveEveryXName = "partialBergmEstimate.rda",
                  seed = NULL,
                  adapt = TRUE,
                  target.accept = 0.25,
                  delayed.acceptance = TRUE,
                  da.inflate = 1.5,
                  cores = 1,
                  log.prior = NULL,
                  prior.decay = c(meanlog = log(0.5), sdlog = 0.75),
                  verbose = TRUE,
                  ...) {

  clock.start <- Sys.time()
  this.call <- match.call()
  .with_seed(seed)

  setup <- .ergm_setup(formula, constraints, control.args = list(...))
  p <- setup$p
  if (network.naedgecount(setup$network) > 0)
    message("Network has missing data. Use bergmM() instead.")
  .check_offsets(setup, offset.coef)

  mp <- .mple(formula, setup, offset.coef, constraints)
  mple <- mp$coef
  fx <- .fix_infinite(mple, setup$offset, offset.coef, verbose)
  offset <- fx$offset
  offset.coef <- fx$offset.coef
  free <- !offset
  pf <- sum(free)
  if (pf < 1L) stop("The model has no free parameters.", call. = FALSE)

  dec <- .decay_index(setup, offset)
  pm <- .param_map(setup, dec)
  prior <- .make_prior(prior.mean, prior.sigma, p, free, setup$coef.names, log.prior,
                       dec = dec, prior.decay = prior.decay)
  prev <- NULL
  if (inherits(startVals, "bergm")) {
    prev <- startVals
    startVals <- NULL
    if (is.null(prev$chains) || !identical(unname(prev$specs), unname(setup$coef.names)))
      stop("'startVals' is a bergm fit of a different model or without per-chain output.",
           call. = FALSE)
    if (is.null(nchains)) nchains <- length(prev$chains)
    if (missing(gamma) && !is.null(prev$gamma)) gamma <- prev$gamma
  }
  blocked <- cores > 1
  if (is.null(nchains)) {
    nchains <- max(2L * pf, 4L)
    if (blocked) nchains <- ceiling(nchains / (2 * cores)) * 2 * cores
  }
  nchains <- as.integer(nchains)
  if (nchains < 3L) stop("'nchains' must be at least 3.", call. = FALSE)
  if (blocked && nchains < 4L) stop("Parallel ADS needs at least 4 chains.", call. = FALSE)

  if (!is.null(prev)) {
    last <- vapply(prev$chains, function(ch) as.matrix(ch)[nrow(as.matrix(ch)), ], numeric(p))
    theta0 <- matrix(last, nrow = p)[, rep_len(seq_len(ncol(matrix(last, nrow = p))), nchains),
                                     drop = FALSE]
    dimnames(theta0) <- list(setup$coef.names, NULL)
    theta0 <- pm$to_work(theta0)
  } else {
    theta0 <- .initial_population(p, nchains, free, offset.coef, mple, mp$vcov,
                                  startVals, setup$coef.names, offset,
                                  .valid_logprior(prior$logdens, pm$eta),
                                  dec = dec, decay = prior$decay)
  }
  U_noise <- .chol_pd(.as_cov(V.proposal, p, "V.proposal")[free, free, drop = FALSE])

  st <- ergm_state_send(setup$state)
  sim_one <- function(h, eta) {
    z <- .simulate_state(st, eta, aux.iters)
    if (is.null(z)) NULL else z$s[1L, ]
  }
  sim_many <- NULL
  if (blocked) {
    cl <- .start_cluster(cores, seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, .worker_set, "state", st)
    sim_many <- function(hs, etas) {
      chunks <- .chunk(seq_along(etas), length(cl))
      res <- parallel::clusterApply(cl, lapply(chunks, function(ix) etas[ix]),
                                    .worker_sim, burnin = aux.iters)
      lapply(unlist(res, recursive = FALSE), function(m) if (is.null(m)) NULL else m[1L, ])
    }
  }

  checkpoint <- NULL
  if (!is.null(saveEveryX))
    checkpoint <- list(at = seq(saveEveryX, burn.in + main.iters, by = saveEveryX),
                       file = saveEveryXName, extra = function() list())

  if (verbose) message(" > MCMC start")
  res <- .exchange_engine(theta = theta0, free = free, eta_fn = pm$eta, out_fn = pm$to_nat,
                          logprior = prior$logdens,
                          sy = matrix(setup$sy, setup$q, nchains),
                          sim_one = sim_one, sim_many = sim_many,
                          burn.in = burn.in, main.iters = main.iters, thin = thin,
                          gamma = gamma, U_noise = U_noise, adapt = adapt,
                          target.accept = target.accept,
                          da = if (isTRUE(delayed.acceptance) &&
                                   (burn.in >= 10 || !missing(delayed.acceptance)))
                            list(inflate = da.inflate),
                          checkpoint = checkpoint, blocked = blocked, verbose = verbose)
  runtime <- difftime(Sys.time(), clock.start)

  .bergm_output(res, setup, offset, prior, main.iters, thin, cut.reject, formula, runtime,
                extra = list(mple = mple, call = this.call,
                             settings = list(burn.in = burn.in, main.iters = main.iters,
                                             aux.iters = aux.iters, thin = thin,
                                             cores = cores, adapt = adapt,
                                             delayed.acceptance = delayed.acceptance)))
}

.valid_logprior <- function(logprior, eta_fn) {
  function(th) if (anyNA(eta_fn(th))) -Inf else logprior(th)
}

.initial_population <- function(p, nchains, free, offset.coef, mple, V, startVals,
                                coef.names, offset, logprior = NULL,
                                dec = rep(FALSE, p), decay = NULL) {
  draw <- function(n, spread) {
    th <- .draw_population(p, n, free, offset.coef, mple, V, startVals, coef.names,
                           offset, spread, dec)
    if (any(dec)) {
      if (!is.null(startVals)) {
        th[dec, ] <- log(startVals[dec]) + stats::runif(sum(dec) * n, -0.1, 0.1) * spread
      } else {
        b <- mple[dec]
        b <- ifelse(is.finite(b) & b > 0.05 & b < 5, log(b), decay$meanlog)
        th[dec, ] <- b + stats::rnorm(sum(dec) * n, 0, pmin(0.5, decay$sdlog)) * spread
      }
    }
    th
  }
  if (any(dec) && !is.null(startVals) && any(startVals[dec] <= 0))
    stop("Decay parameters in 'startVals' must be positive.", call. = FALSE)
  theta <- draw(nchains, 1)
  if (is.null(logprior)) return(theta)
  for (attempt in seq_len(200)) {
    bad <- which(!is.finite(apply(theta, 2, logprior)))
    if (!length(bad)) return(theta)
    theta[, bad] <- draw(length(bad), 1 + attempt / 20)
  }
  stop("Could not find starting values with positive prior density; supply 'startVals'.",
       call. = FALSE)
}

.draw_population <- function(p, nchains, free, offset.coef, mple, V, startVals,
                             coef.names, offset, spread = 1, dec = rep(FALSE, p)) {
  theta <- matrix(NA_real_, p, nchains, dimnames = list(coef.names, NULL))
  g <- free & !dec
  if (!is.null(startVals)) {
    if (length(startVals) != p)
      stop(sprintf("'startVals' must have length %d.", p), call. = FALSE)
    theta[] <- startVals + stats::runif(p * nchains, -0.1, 0.1) * spread
  } else if (any(g)) {
    base <- mple
    base[!is.finite(base)] <- 0
    Vf <- if (!is.null(V) && all(dim(V) == p)) V[g, g, drop = FALSE] else NULL
    if (!is.null(Vf) && all(is.finite(Vf)) && .is_pd(Vf)) {
      f <- pmin(1, 1 / sqrt(diag(Vf)))
      Vf <- Vf * outer(f, f) * spread^2
      theta[g, ] <- t(.rmvn(nchains, base[g], chol(Vf)))
    } else {
      theta[g, ] <- base[g] + stats::runif(sum(g) * nchains, -0.1, 0.1) * spread
    }
  }
  if (any(offset)) theta[offset, ] <- offset.coef
  theta
}
