#' Adjustment of ERGM pseudolikelihood
#'
#' Function to estimate the transformation parameters for
#' adjusting the pseudolikelihood function (Bouranis, Friel and Maire, 2018).
#'
#' The log-likelihood at the MLE, and hence the adjustment constant, is
#' estimated by default (\code{llk.method = "bridge"}) by bridge sampling from
#' the dyad-independent submodel of the formula, using
#' \code{\link[ergm]{ergm.bridge.dindstart.llk}}. The reference model has an
#' exact likelihood and is close to the target, which keeps the estimate
#' accurate on large networks. \code{llk.method = "ladder"} uses the
#' stepping-stone ladder \eqn{\theta_t = t\,\hat\theta_{MLE}}, \eqn{t \in [0,1]},
#' of previous versions, traversed in both directions with warm-started chains
#' and averaging forward and reverse estimates (which removes most of the bias
#' of a single sweep); it needs many rungs when the network is large, because
#' the path starts at the uniform random graph. With \code{cores > 1} the two
#' sweeps run in parallel. The returned object can be passed to
#' \code{\link{evidenceCJ}} and \code{\link{evidencePP}} through their
#' \code{apl} argument to avoid recomputing it.
#'
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}.
#'
#' @param aux.iters count; number of auxiliary iterations used for drawing the first network from the ERGM likelihood
#' (for the curvature matching and, with \code{llk.method = "ladder"}, at every rung of the ladder). Default 3000.
#'
#' @param n.aux.draws count; number of auxiliary networks drawn from the ERGM likelihood. Default 50.
#'
#' @param aux.thin count; number of auxiliary iterations between network draws after the first network is drawn. Default 50.
#'
#' @param ladder count; length of temperature ladder (>=3) for \code{llk.method = "ladder"}. Default 50.
#'
#' @param estimate If "MLE" (the default), then an approximate maximum likelihood estimator is returned.
#' If "CD", the Monte-Carlo contrastive divergence estimate is returned. See \code{\link[ergm]{ergm}}.
#'
#' @param seed integer; seed for the random number generator. See \code{set.seed}.
#'
#' @param mle optional; a fitted \code{\link[ergm]{ergm}} object or a numeric vector of
#' coefficients to use as the MLE, skipping the \code{ergm()} fit.
#'
#' @param cores count; number of parallel workers for the ladder (at most 2 are used,
#' one per sweep).
#'
#' @param llk.method character; \code{"bridge"} (default) or \code{"ladder"}, the method
#' used to estimate the log-likelihood at the MLE (see Details).
#'
#' @param bridge.control optional \code{\link[ergm]{control.ergm.bridge}} list for
#' \code{llk.method = "bridge"}.
#'
#' @param ... Additional arguments, to be passed to the ergm function. See \code{\link[ergm]{ergm}}.
#'
#' @return An object of class \code{ergmAPL}: a list with \code{formula},
#' \code{Theta_MLE}, \code{Theta_PL}, \code{W}, \code{logC}, \code{ll_true},
#' \code{logztheta} and \code{ll_adjpseudo}.
#'
#' @references
#' Bouranis, L., Friel, N., & Maire, F. (2018). Bayesian model selection for exponential
#' random graph models via adjusted pseudolikelihoods.
#' Journal of Computational and Graphical Statistics, 27(3), 516-528. \url{https://arxiv.org/abs/1706.06344}
#'
#' @examples
#' \dontrun{
#' data(florentine)
#' apl <- ergmAPL(flomarriage ~ edges + kstar(2))
#' apl
#' }
#' @export
ergmAPL <- function(formula,
                    aux.iters   = NULL,
                    n.aux.draws = NULL,
                    aux.thin    = NULL,
                    ladder      = NULL,
                    estimate    = c("MLE", "CD"),
                    seed        = 1,
                    mle         = NULL,
                    cores       = 1,
                    llk.method  = c("bridge", "ladder"),
                    bridge.control = NULL,
                    ...) {
  clock.start <- Sys.time()
  estimate <- match.arg(estimate)
  llk.method <- match.arg(llk.method)
  aux.iters <- aux.iters %||% 3000
  n.aux.draws <- n.aux.draws %||% 50
  aux.thin <- aux.thin %||% 50
  ladder <- ladder %||% 50
  if (ladder < 3) stop("'ladder' must be at least 3.", call. = FALSE)
  .with_seed(seed)

  setup <- .ergm_setup(formula)
  .not_supported_curved(setup, "ergmAPL")
  y <- setup$network
  sy <- setup$sy
  specs <- setup$coef.names
  pl <- .pl_setup(formula)

  theta_mle <- if (is.null(mle)) {
    .ergm_estimate(formula, estimate, seed, ...)
  } else if (inherits(mle, "ergm")) stats::coef(mle) else as.numeric(mle)
  names(theta_mle) <- specs
  theta_pl <- .pl_mode(pl, init = theta_mle)
  if (is.null(theta_pl)) {
    theta_pl <- tryCatch(stats::coef(suppressMessages(ergm(formula, estimate = "MPLE"))),
                         error = function(e) rep(Inf, length(specs)))
  }
  names(theta_pl) <- specs

  if (any(!is.finite(theta_mle)) || any(!is.finite(theta_pl))) {
    inf_term <- specs[!is.finite(theta_mle) | !is.finite(theta_pl)]
    stop("The MLE/MPLE contains an infinite value for the following model terms: ",
         paste(inf_term, collapse = " & "), ". Consider changing these model terms.",
         call. = FALSE)
  }

  st <- ergm_state_send(setup$state)
  if (isTRUE(tryCatch(is.dyad.independent(setup$model), error = function(e) FALSE))) {
    W <- diag(length(specs))
    dimnames(W) <- list(specs, specs)
  } else {
    z <- .simulate_state(st, theta_mle, burnin = aux.iters, samplesize = n.aux.draws,
                         interval = aux.thin)
    if (is.null(z)) stop(.sim_limit_msg, call. = FALSE)
    W <- .apl_W(stats::cov(z$s), .pl_hess(pl, theta_pl))
  }
  adjust.info <- list(Theta_MLE = theta_mle, Theta_PL = theta_pl, W = W)

  logz0A <- network.dyadcount(y, na.omit = FALSE) * log(2)
  logE <- NULL
  if (llk.method == "bridge") {
    bc <- bridge.control %||% control.ergm.bridge(seed = seed)
    ll.true <- suppressMessages(suppressWarnings(
      ergm.bridge.dindstart.llk(formula, coef = theta_mle, control = bc)))
    ll.true <- as.numeric(ll.true)
    logztheta <- sum(theta_mle * sy) - ll.true
  } else {
    path <- seq(0, 1, length.out = ladder)
    sweeps <- if (cores > 1) {
      cl <- .start_cluster(min(cores, 2L), seed)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterApply(cl, c("up", "down"), .apl_sweep, st = st, path = path,
                             theta = theta_mle, burnin = aux.iters, n = n.aux.draws,
                             interval = aux.thin)
    } else {
      lapply(c("up", "down"), .apl_sweep, st = st, path = path, theta = theta_mle,
             burnin = aux.iters, n = n.aux.draws, interval = aux.thin)
    }
    logE <- .apl_ladder_logratios(sweeps, path)
    logztheta <- sum(logE) + logz0A
    ll.true <- sum(theta_mle * sy) - logztheta
  }
  ll.adjpseudo <- .pl_loglik(pl, theta_pl)

  out <- list(formula      = formula,
              Theta_MLE    = theta_mle,
              Theta_PL     = theta_pl,
              W            = W,
              logC         = ll.true - ll.adjpseudo,
              ll_true      = ll.true,
              logztheta    = logztheta,
              ll_adjpseudo = ll.adjpseudo,
              ladder.logE  = logE,
              llk.method   = llk.method,
              Time         = difftime(Sys.time(), clock.start))
  class(out) <- "ergmAPL"
  out
}

.apl_sweep <- function(direction, st, path, theta, burnin, n, interval) {
  idx <- if (direction == "up") seq_along(path) else rev(seq_along(path))
  U <- vector("list", length(path))
  state <- st
  for (i in idx) {
    z <- .simulate_state(state, path[i] * theta, burnin = burnin, samplesize = n,
                         interval = interval)
    if (is.null(z)) stop(.sim_limit_msg, call. = FALSE)
    state <- update(state, state = z$state)
    U[[i]] <- drop(z$s %*% theta)
  }
  U
}

.apl_ladder_logratios <- function(sweeps, path) {
  dt <- diff(path)
  vapply(seq_along(dt), function(i) {
    est <- unlist(lapply(sweeps, function(U)
      c(.logmeanexp(dt[i] * U[[i]]), -.logmeanexp(-dt[i] * U[[i + 1]]))))
    mean(est)
  }, numeric(1))
}

.apl_W <- function(covS, H_pl) {
  A <- .make_pd(covS)
  B <- .make_pd(-H_pl)
  W <- backsolve(chol(B), chol(A))
  dimnames(W) <- dimnames(H_pl)
  W
}

#' @export
print.ergmAPL <- function(x, ...) {
  cat("Adjusted pseudolikelihood for model: y ~", .formula_rhs(x$formula), "\n\n")
  print(rbind(MLE = x$Theta_MLE, MPLE = x$Theta_PL))
  cat("\nlog normalising constant at the MLE:", format(x$logztheta),
      "\nlog-likelihood at the MLE:", format(x$ll_true),
      "\nlog C (adjustment constant):", format(x$logC), "\n")
  invisible(x)
}
