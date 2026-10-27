#' Evidence estimation via Chib and Jeliazkov's method
#'
#' Function to estimate the evidence (marginal likelihood) with Chib and Jeliazkov's method,
#' based on the adjusted pseudolikelihood function.
#'
#' The adjusted pseudo-posterior is sampled with a compiled random-walk
#' Metropolis kernel and all pseudolikelihood evaluations needed by the
#' estimator are vectorised in compiled code. The Metropolis-Hastings
#' acceptance probabilities in the Chib-Jeliazkov identity are computed with
#' the full (pseudo-)posterior, i.e. including the prior.
#'
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}.
#'
#' @param prior.mean vector;
#' mean vector of the multivariate Normal prior.
#' By default set to a vector of 0's.
#'
#' @param prior.sigma square matrix (or scalar/vector diagonal);
#' variance/covariance matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100.
#'
#' @param aux.iters count;
#' number of auxiliary iterations used for drawing the first network from the ERGM likelihood.
#' See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param n.aux.draws count;
#' number of auxiliary networks drawn from the ERGM likelihood.
#' See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param aux.thin count;
#' number of auxiliary iterations between network draws after the first network is drawn.
#' See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param ladder count; length of temperature ladder (>=3).
#' See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param main.iters count;
#' number of MCMC iterations after burn-in for the adjusted pseudo-posterior estimation.
#'
#' @param burn.in count;
#' number of burn-in iterations at the beginning of an MCMC run
#' for the adjusted pseudo-posterior estimation.
#'
#' @param thin count;
#' thinning interval used in the simulation for the adjusted pseudo-posterior estimation.
#'
#' @param num.samples integer;
#' number of proposal draws used in the marginal likelihood estimate.
#'
#' @param V.proposal count;
#' diagonal entry for the multivariate Normal proposal.
#' By default set to 1.5.
#'
#' @param seed integer; seed for the random number generator.
#'
#' @param estimate If "MLE" (the default), then an approximate maximum likelihood estimator is returned.
#' If "CD" , the Monte-Carlo contrastive divergence estimate is returned. See \code{\link[ergm]{ergm}}.
#'
#' @param apl optional; an \code{ergmAPL} object for the same model, to reuse
#' an existing pseudolikelihood adjustment.
#'
#' @param cores count; number of parallel workers used by \code{\link{ergmAPL}}.
#'
#' @param verbose logical; print progress information.
#'
#' @param ... additional arguments, to be passed to the ergm function.
#' See \code{\link[ergm]{ergm}} and \code{\link[Bergm]{ergmAPL}}.
#'
#' @references
#' Caimo, A., & Friel, N. (2013). Bayesian model selection for exponential random graph models.
#' Social Networks, 35(1), 11-24. \url{https://arxiv.org/abs/1201.2337}
#'
#' Bouranis, L., Friel, N., & Maire, F. (2018). Bayesian model selection for exponential
#' random graph models via adjusted pseudolikelihoods.
#' Journal of Computational and Graphical Statistics, 27(3), 516-528.
#' \url{https://arxiv.org/abs/1706.06344}
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network:
#' data(florentine)
#'
#' # MCMC sampling and evidence estimation:
#' CJE <- evidenceCJ(flomarriage ~ edges + kstar(2),
#'                   main.iters  = 2000,
#'                   burn.in     = 200,
#'                   aux.iters   = 500,
#'                   num.samples = 25000,
#'                   V.proposal  = 2.5)
#'
#' # Posterior summaries:
#' summary(CJE)
#'
#' # MCMC diagnostics plots:
#' plot(CJE)
#'
#' # Log-evidence (marginal likelihood) estimate:
#' CJE$log.evidence
#'}
#'
#' @export
evidenceCJ <- function(formula,
                       prior.mean  = NULL,
                       prior.sigma = NULL,
                       aux.iters   = 1000,
                       n.aux.draws = 50,
                       aux.thin    = 50,
                       ladder      = 30,
                       main.iters  = 30000,
                       burn.in     = 5000,
                       thin        = 1,
                       V.proposal  = 1.5,
                       num.samples = 25000,
                       seed        = 1,
                       estimate    = c("MLE", "CD"),
                       apl         = NULL,
                       cores       = 1,
                       verbose     = TRUE,
                       ...) {
  clock.start <- Sys.time()
  estimate <- match.arg(estimate)

  if (is.null(apl)) {
    if (verbose) message(" > Pseudolikelihood adjustment")
    apl <- ergmAPL(formula = formula, aux.iters = aux.iters, n.aux.draws = n.aux.draws,
                   aux.thin = aux.thin, ladder = ladder, estimate = estimate,
                   seed = seed, cores = cores, ...)
  }
  .with_seed(seed)
  pl <- .pl_setup(formula)
  dim <- pl$p
  specs <- names(apl$Theta_MLE) %||% pl$names
  prior <- .gauss_prior(prior.mean, prior.sigma, dim, specs)

  H <- t(apl$W) %*% .pl_hess(pl, apl$Theta_PL) %*% apl$W
  Sigma.proposal <- .as_cov(V.proposal, dim, "V.proposal")
  S.prop <- Sigma.proposal %*% solve(prior$prec - H) %*% Sigma.proposal
  S.prop <- (S.prop + t(S.prop)) / 2

  if (verbose) message(" > MCMC run")
  run <- .rwm_pl(pl, prior, S.prop, init = apl$Theta_MLE, burn.in = burn.in,
                 main.iters = main.iters, thin = thin, adj = apl)
  T0 <- coda::mcmc(run$draws, thin = thin)

  if (verbose) message(" > Model evidence estimation")
  logpost <- function(Theta) {
    if (is.null(dim(Theta))) Theta <- matrix(Theta, nrow = 1)
    .adjpl_batch(pl, Theta, apl)$logpl + prior$logdens(Theta)
  }
  thetahat <- colMeans(run$draws)
  log.post <- logpost(thetahat)

  draws <- run$draws
  theta.g <- if (num.samples >= nrow(draws)) draws else
    draws[sample.int(nrow(draws), num.samples), , drop = FALSE]
  U.prop <- .chol_pd(S.prop)
  q.g <- exp(.dmvn_log(theta.g, thetahat, U.prop))
  alpha.g <- pmin(1, exp(log.post - logpost(theta.g)))
  theta.j <- .rmvn(num.samples, thetahat, U.prop)
  alpha.j <- pmin(1, exp(logpost(theta.j) - log.post))
  pi.hat <- mean(alpha.g * q.g) / mean(alpha.j)
  logEvidence <- apl$logC + log.post - log(pi.hat)

  ess <- round(coda::effectiveSize(T0), 0)
  names(ess) <- specs
  out <- list(formula      = formula,
              Theta        = T0,
              AR           = round(run$accept, 2),
              ess          = ess,
              log.evidence = unname(logEvidence),
              dim          = dim,
              specs        = specs,
              fixed        = stats::setNames(rep(FALSE, dim), specs),
              Time         = difftime(Sys.time(), clock.start),
              apl          = apl)
  class(out) <- "bergm"
  out
}
