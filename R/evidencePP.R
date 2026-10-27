#' Evidence estimation via power posteriors
#'
#' Function to estimate the evidence (marginal likelihood) with Power posteriors,
#' based on the adjusted pseudolikelihood function.
#'
#' Each tempered adjusted pseudo-posterior is sampled with a compiled
#' random-walk Metropolis kernel, and the zero-variance control variates and
#' the corrected trapezoidal rule (Friel, Hurn and Wyse, 2014) are computed
#' with vectorised compiled code.
#'
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}.
#'
#' @param prior.mean vector; mean vector of the multivariate Normal prior.
#' By default set to a vector of 0's.
#'
#' @param prior.sigma square matrix (or scalar/vector diagonal); variance/covariance matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100.
#'
#' @param aux.iters count; number of auxiliary iterations used for drawing the first network from the ERGM likelihood.
#' See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param n.aux.draws count; number of auxiliary networks drawn from the ERGM likelihood.
#' See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param aux.thin count; number of auxiliary iterations between network draws after the first network is drawn.
#' See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param ladder count; length of temperature ladder (>=3). See \code{\link[Bergm]{ergmAPL}}.
#'
#' @param main.iters count; number of MCMC iterations after burn-in for the adjusted pseudo-posterior estimation.
#'
#' @param burn.in count; number of burn-in iterations at the beginning of an MCMC run for the adjusted pseudo-posterior estimation.
#'
#' @param thin count; thinning interval used in the simulation for the adjusted pseudo-posterior estimation.
#'
#' @param V.proposal count; diagonal entry for the multivariate Normal proposal.
#' By default set to 1.5.
#'
#' @param seed integer; seed for the random number generator.
#'
#' @param temps numeric vector; inverse temperature ladder, \eqn{t \in [0,1]}.
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
#' Bouranis, L., Friel, N., & Maire, F. (2018). Bayesian model selection for exponential
#' random graph models via adjusted pseudolikelihoods.
#' Journal of Computational and Graphical Statistics, 27(3), 516-528.
#' \url{https://arxiv.org/abs/1706.06344}
#'
#' Friel, N., Hurn, M., & Wyse, J. (2014). Improving power posterior estimation of
#' statistical evidence. Statistics and Computing, 24(5), 709-723.
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network:
#' data(florentine)
#'
#' PPE <- evidencePP(flomarriage ~ edges + kstar(2),
#'                   aux.iters   = 500,
#'                   aux.thin    = 50,
#'                   main.iters  = 2000,
#'                   burn.in     = 100,
#'                   V.proposal  = 2.5)
#'
#' # Posterior summaries:
#' summary(PPE)
#'
#' # MCMC diagnostics plots:
#' plot(PPE)
#'
#' # Log-evidence (marginal likelihood) estimate:
#' PPE$log.evidence
#'}
#'
#' @export
evidencePP <- function(formula,
                       prior.mean  = NULL,
                       prior.sigma = NULL,
                       aux.iters   = 1000,
                       n.aux.draws = 50,
                       aux.thin    = 50,
                       ladder      = 30,
                       main.iters  = 20000,
                       burn.in     = 5000,
                       thin        = 1,
                       V.proposal  = 1.5,
                       seed        = 1,
                       temps       = NULL,
                       estimate    = c("MLE", "CD"),
                       apl         = NULL,
                       cores       = 1,
                       verbose     = TRUE,
                       ...) {
  clock.start <- Sys.time()
  estimate <- match.arg(estimate)
  if (is.null(temps)) temps <- seq(0, 1, length.out = 50)^5
  temps <- sort(as.numeric(temps))
  numtemp <- length(temps)
  if (numtemp < 3 || temps[1] != 0 || temps[numtemp] != 1)
    stop("'temps' must be an increasing ladder from 0 to 1 with at least 3 values.",
         call. = FALSE)

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
  dimnames(S.prop) <- list(specs, specs)
  taup <- 1 / sqrt(diag(S.prop))
  tau0 <- prior$prec
  a.tilde <- log(tau0[1, 1] / taup) / log(temps[2])
  covar.temps <- vector("list", numtemp)
  for (i in numtemp:1) {
    covar.temps[[i]] <- S.prop
    diag(covar.temps[[i]]) <- diag(S.prop) / (temps[i]^a.tilde[1])
  }
  covar.temps[[1]] <- covar.temps[[2]]

  if (verbose) message(" > MCMC run")
  pplist <- vector("list", numtemp)
  acceptances <- numeric(numtemp)
  l0 <- apl$Theta_MLE
  for (l in numtemp:1) {
    if (l == 1) {
      pplist[[l]] <- .rmvn(main.iters, prior$mean, prior$U)
      colnames(pplist[[l]]) <- specs
      acceptances[l] <- 1
    } else {
      run <- .rwm_pl(pl, prior, .make_pd(covar.temps[[l]]), init = l0, burn.in = burn.in,
                     main.iters = main.iters, thin = thin, temperature = temps[l], adj = apl)
      pplist[[l]] <- run$draws
      acceptances[l] <- round(run$accept, 3)
      l0 <- run$last
    }
  }

  if (verbose) message(" > Model evidence estimation")
  names(acceptances) <- paste("Temp", "_", 1:numtemp, " = ", round(temps, 5), sep = "")
  cv.ell.D2 <- cv.vll.D2 <- numeric(numtemp)
  for (j in seq_len(numtemp)) {
    Th <- pplist[[j]]
    ev <- .adjpl_batch(pl, Th, apl, grad = TRUE)
    prior_score <- -t(prior$prec %*% (t(Th) - prior$mean))
    grad_lpp <- temps[j] * ev$grad + prior_score
    cv <- .zv_control_variate(Th, grad_lpp, ev$logpl)
    cv.ell.D2[j] <- cv$mean
    cv.vll.D2[j] <- cv$var
  }

  pp.estimates <- .ppml(cv.ell.D2, cv.vll.D2, temps)
  log.evidence <- apl$logC + pp.estimates
  Theta <- coda::mcmc(pplist[[numtemp]], thin = thin)
  ess <- round(coda::effectiveSize(Theta), 0)
  names(ess) <- specs

  out <- list(formula      = formula,
              Theta        = Theta,
              AR_temp      = acceptances,
              AR           = acceptances[[numtemp]],
              ess          = ess,
              log.evidence = unname(log.evidence),
              dim          = dim,
              specs        = specs,
              fixed        = stats::setNames(rep(FALSE, dim), specs),
              Time         = difftime(Sys.time(), clock.start),
              apl          = apl)
  class(out) <- "bergm"
  out
}

.zv_control_variate <- function(Th, g, f) {
  d <- ncol(Th)
  cols <- list(g, Th * g + 1)
  if (d > 1) {
    pr <- utils::combn(d, 2)
    cols[[3]] <- vapply(seq_len(ncol(pr)), function(k) {
      a <- pr[1, k]; b <- pr[2, k]
      Th[, a] * g[, b] + Th[, b] * g[, a]
    }, numeric(nrow(Th)))
  }
  w.mat <- do.call(cbind, cols)
  V <- stats::var(w.mat)
  cvec <- stats::cov(w.mat, f)
  phi <- tryCatch(-solve(V, cvec), error = function(e) {
    co <- -qr.coef(qr(V), cvec)
    co[is.na(co)] <- 0
    co
  })
  h <- f + drop(w.mat %*% phi)
  list(mean = mean(h), var = stats::var(h))
}

.ppml <- function(cv.ell.D2, cv.vll.D2, tempvec) {
  N <- length(cv.ell.D2)
  res <- 0
  for (i in 1:(N - 1)) {
    wts <- tempvec[i + 1] - tempvec[i]
    res <- res + wts * ((cv.ell.D2[i + 1] + cv.ell.D2[i]) / 2) -
      ((wts^2) / 12) * (cv.vll.D2[i + 1] - cv.vll.D2[i])
  }
  res
}
