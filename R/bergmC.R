#' Calibrating misspecified Bayesian ERGMs
#'
#' Function to transform a sample from the pseudo-posterior
#' to one that is approximately sampled from the intractable
#' posterior distribution (Bouranis, Friel and Maire, 2017).
#'
#' The pseudo-posterior is sampled by a compiled random-walk Metropolis kernel.
#' The mode of the true posterior, which anchors the calibration, is found by
#' default with a Newton-type stochastic approximation that uses the
#' simulated covariance of the sufficient statistics and a persistent
#' auxiliary network chain (\code{map.method = "newton"}); this typically
#' converges in about ten iterations, whereas the Robbins-Monro recursion of
#' previous versions (\code{map.method = "rm"}) needs hundreds. The
#' pseudo-posterior mode is found exactly by Newton-Raphson.
#'
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}.
#'
#' @param prior.mean vector; mean vector of the multivariate Normal prior.
#' By default set to a vector of 0's.
#'
#' @param prior.sigma square matrix (or scalar/vector diagonal); variance/covariance
#' matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100.
#'
#' @param burn.in count; number of burn-in iterations at the beginning of an MCMC run for the pseudo-posterior estimation.
#'
#' @param main.iters count; number of MCMC iterations after burn-in for the pseudo-posterior estimation.
#'
#' @param aux.iters count; number of auxiliary iterations used for drawing the first network from the ERGM likelihood.
#'
#' @param V.proposal count; diagonal entry for the multivariate Normal proposal.
#' By default set to 1.5.
#'
#' @param thin count; thinning interval used in the simulation for the pseudo-posterior estimation.
#'
#' @param rm.iters count; number of iterations for the Robbins-Monro algorithm (\code{map.method = "rm"}).
#'
#' @param rm.a scalar; constant for sequence alpha_n (Robbins-Monro).
#'
#' @param rm.alpha scalar; noise added to gradient (Robbins-Monro).
#'
#' @param n.aux.draws count; number of auxiliary networks drawn from the ERGM likelihood at every
#' stochastic approximation iteration and for the curvature adjustment.
#'
#' @param aux.thin count; number of auxiliary iterations between network draws after the first network is drawn.
#'
#' @param estimate If "MLE" (the default), then an approximate maximum likelihood estimator is used as a starting point
#' for the posterior mode search. If "CD", the Monte-Carlo contrastive divergence estimate is used. See \code{\link[ergm]{ergm}}.
#'
#' @param seed integer; seed for the random number generator. See \code{set.seed}.
#'
#' @param map.method character; \code{"newton"} (default) or \code{"rm"} (Robbins-Monro, as in Bergm 5).
#'
#' @param map.iters count; number of Newton stochastic approximation iterations.
#'
#' @param init vector; optional starting value for the posterior mode search. If given,
#' the \code{\link[ergm]{ergm}} fit is skipped.
#'
#' @param verbose logical; print progress information.
#'
#' @param ... Additional arguments, to be passed to the ergm function. See \code{\link[ergm]{ergm}}.
#'
#' @return An object of class \code{bergm}, with additional elements
#' \code{Theta_star} (estimated posterior mode), \code{Theta_PL}
#' (pseudo-posterior mode), \code{W} (curvature adjustment matrix) and
#' \code{map.trajectory}.
#'
#' @references
#' Bouranis, L., Friel, N., & Maire, F. (2017). Efficient Bayesian inference for exponential
#' random graph models by correcting the pseudo-posterior distribution.
#' Social Networks, 50, 98-108. \url{https://arxiv.org/abs/1510.00934}
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data(florentine)
#'
#' # Calibrated pseudo-posterior:
#' cpp.flo <- bergmC(flomarriage ~ edges + kstar(2),
#'                   aux.iters  = 500,
#'                   burn.in    = 500,
#'                   main.iters = 10000,
#'                   V.proposal = 2.5)
#'
#' # Posterior summaries:
#' summary(cpp.flo)
#'}
#'
#' @export
bergmC <- function(formula,
                   prior.mean  = NULL,
                   prior.sigma = NULL,
                   burn.in     = 1e04,
                   main.iters  = 4e04,
                   aux.iters   = 3000,
                   V.proposal  = 1.5,
                   thin        = 1,
                   rm.iters    = 500,
                   rm.a        = 0.001,
                   rm.alpha    = 0,
                   n.aux.draws = 400,
                   aux.thin    = 50,
                   estimate    = c("MLE", "CD"),
                   seed        = 1,
                   map.method  = c("newton", "rm"),
                   map.iters   = 20,
                   init        = NULL,
                   verbose     = TRUE,
                   ...) {
  clock.start <- Sys.time()
  estimate <- match.arg(estimate)
  map.method <- match.arg(map.method)
  .with_seed(seed)

  setup <- .ergm_setup(formula)
  .not_supported_curved(setup, "bergmC")
  dim <- setup$p
  specs <- setup$coef.names
  sy <- setup$sy
  prior <- .gauss_prior(prior.mean, prior.sigma, dim, specs)
  pl <- .pl_setup(formula)

  mode_pl <- .pl_mode(pl, prior = prior)
  if (is.null(mode_pl)) stop("Could not locate the pseudo-posterior mode.", call. = FALSE)
  H_pl <- .pl_hess(pl, mode_pl)
  Sigma.proposal <- .as_cov(V.proposal, dim, "V.proposal")
  S.prop <- Sigma.proposal %*% solve(prior$prec - H_pl) %*% Sigma.proposal

  if (verbose) message(" > MCMC start")
  unadj <- .rwm_pl(pl, prior, S.prop, init = mode_pl, burn.in = burn.in,
                   main.iters = main.iters, thin = thin)

  if (verbose) message(" > MAP estimation")
  if (is.null(init)) init <- .ergm_estimate(formula, estimate, seed, ...)
  init[!is.finite(init)] <- 0
  st <- ergm_state_send(setup$state)
  map <- if (map.method == "newton") {
    .map_newton(st, sy, init, prior, aux.iters, n.aux.draws, aux.thin, iters = map.iters)
  } else {
    .map_rm(st, sy, init, prior, aux.iters, n.aux.draws, aux.thin, rm.iters, rm.a, rm.alpha)
  }
  theta_star <- stats::setNames(map$theta, specs)
  if (map.method == "newton" && isTRUE(map$discrepancy > 2))
    warning("The posterior mode search did not match the observed statistics well ",
            "(Mahalanobis discrepancy ", round(map$discrepancy, 2), "). The model may be ",
            "near-degenerate, or 'aux.iters'/'n.aux.draws' too small.", call. = FALSE)

  if (verbose) message(" > Curvature Adjustment")
  z <- .simulate_state(map$state, theta_star, burnin = aux.iters,
                       samplesize = n.aux.draws, interval = aux.thin)
  if (is.null(z)) stop(.sim_limit_msg, call. = FALSE)
  chol.true.posterior <- .chol_pd(stats::cov(z$s) + prior$prec)
  chol.PL.posterior <- .chol_pd(prior$prec - H_pl)
  W <- backsolve(chol.PL.posterior, chol.true.posterior)
  V <- solve(W)
  corrected <- sweep(sweep(unadj$draws, 2, mode_pl) %*% t(V), 2, theta_star, "+")
  colnames(corrected) <- specs

  runtime <- difftime(Sys.time(), clock.start)
  FF <- coda::mcmc(corrected, thin = thin)
  ess <- round(coda::effectiveSize(FF), 0)
  names(ess) <- specs
  fixed <- stats::setNames(rep(FALSE, dim), specs)

  out <- list(Theta_star = theta_star,
              Theta      = FF,
              Time       = runtime,
              formula    = formula,
              AR         = round(unadj$accept, 2),
              ess        = ess,
              dim        = dim,
              specs      = specs,
              fixed      = fixed,
              Theta_PL   = mode_pl,
              W          = W,
              map.trajectory = map$trajectory,
              prior      = list(mean = prior$mean, sigma = prior$sigma, type = "normal"),
              constraints = setup$constraints)
  class(out) <- "bergm"
  out
}
