#' Evidence estimation via Chib and Jeliazkov's method
#'
#' Function to estimate the evidence (marginal likelihood) with Chib and Jeliazkov's method, 
#' based on the adjusted pseudolikelihood function.
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
#' @param prior.sigma square matrix; 
#' variance/covariance matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100.
#' 
#' @param aux.iters count; 
#' number of auxiliary iterations used for drawing the first network from the ERGM likelihood. 
#' See \code{\link[ergm]{control.simulate.formula}} and \code{\link[Bergm]{ergmAPL}}.
#' 
#' @param n.aux.draws count; 
#' number of auxiliary networks drawn from the ERGM likelihood. 
#' See \code{\link[ergm]{control.simulate.formula}} and \code{\link[Bergm]{ergmAPL}}.
#' 
#' @param aux.thin count; 
#' number of auxiliary iterations between network draws after the first network is drawn. 
#' See \code{\link[ergm]{control.simulate.formula}} and \code{\link[Bergm]{ergmAPL}}.
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
#' The number of MCMC iterations must be divisible by this value.
#' 
#' @param num.samples integer; 
#' number of samples used in the marginal likelihood estimate. 
#' Must be lower than \code{main.iters} - \code{burnin}.
#' 
#' @param V.proposal count; 
#' diagonal entry for the multivariate Normal proposal.
#' By default set to 1.5.
#' 
#' @param seed seed for the random number generator. 
#' See \code{\link[MCMCpack]{MCMCmetrop1R}}. 
#'
#' @param estimate If "MLE" (the default), then an approximate maximum likelihood estimator is returned. If "CD" , the Monte-Carlo contrastive divergence estimate is returned. See \code{\link[ergm]{ergm}}.
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
#' # MCMC diagnostics and posterior summaries:
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
#'
evidenceCJ <- function(formula, 
                       prior.mean  = NULL, 
                       prior.sigma = NULL, 
                       aux.iters   = 1000, 
                       n.aux.draws = 5, 
                       aux.thin    = 50, 
                       ladder      = 30, 
                       main.iters  = 30000, 
                       burn.in     = 5000, 
                       thin        = 1, 
                       V.proposal  = 1.5, 
                       num.samples = 25000, 
                       seed        = NA,
                       estimate    = c("MLE","CD"),
                       ...) 
{
  y     <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  sy    <- summary(formula)
  dim   <- length(sy)
  
  if (dim == 1) 
    stop("Model dimension must be greater than 1")
  if (is.null(prior.mean)) 
    prior.mean <- rep(0, dim)
  if (is.null(prior.sigma)) 
    prior.sigma <- diag(100, dim, dim)
  
  estimate <- match.arg(estimate)
  
  expit <- function(x) exp(x)/(1 + exp(x))
  
  adjusted_logPL <- function(theta, Y, X, weights, calibr.info) {
    theta_transf <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL)
    xtheta       <- c(X %*% theta_transf)
    log.like     <- sum(dbinom(weights * Y, weights, expit(xtheta), log = TRUE))
    return(log.like)
  }
  
  adjusted_logPP <- function(theta, Y, X, weights, prior.mean, 
                             prior.sigma, calibr.info) {
    theta_transf <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL)
    xtheta    <- c(X %*% theta_transf)
    log.like  <- sum(dbinom(weights * Y, weights, expit(xtheta), log = TRUE))
    log.prior <- dmvnorm(theta, prior.mean, prior.sigma, log = TRUE)
    return(log.like + log.prior)
  }
  
  Hessian_adjusted_logPP <- function(theta, X, weights, calibr.info) {
    p <- exp(as.matrix(X) %*% theta)/(1 + exp(as.matrix(X) %*% theta))
    W <- Diagonal(x = as.vector(weights * p * (1 - p)))
    Hessian <- as.matrix(-t(as.matrix(X)) %*% W %*% as.matrix(X))
    return(t(calibr.info$W) %*% Hessian %*% calibr.info$W)
  }
  
  clock.start <- Sys.time()
  
  message(" > Pseudolikelihood adjustment")
  info.adjustPL <- ergmAPL(formula     = formula, 
                           aux.iters   = aux.iters, 
                           n.aux.draws = n.aux.draws, 
                           aux.thin    = aux.thin, 
                           ladder      = ladder,
                           estimate    = estimate, 
                           ...)
  
  calibr.info <- list(Theta_MLE = info.adjustPL$Theta_MLE, 
                      Theta_PL  = info.adjustPL$Theta_PL, 
                      W         = info.adjustPL$W, 
                      logC      = info.adjustPL$logC)
  
  mplesetup <- ergmMPLE(formula)
  data.glm.initial <- cbind(mplesetup$response, 
                            mplesetup$weights, 
                            mplesetup$predictor)
  colnames(data.glm.initial) <- c("responses", "weights", colnames(mplesetup$predictor))
  Sigma.proposal <- diag(rep(V.proposal, dim), dim, dim)
  
  H <- Hessian_adjusted_logPP(calibr.info$Theta_PL, 
                              X           = data.glm.initial[,-c(1, 2)], 
                              weights     = data.glm.initial[, "weights"], 
                              calibr.info = calibr.info)
  S.prop <- Sigma.proposal %*% solve(solve(prior.sigma) - H) %*% Sigma.proposal
  
  message(" > MCMC run")
  capture.output(T0 <- MCMCmetrop1R(adjusted_logPP, 
                                    theta.init  = calibr.info$Theta_PL, 
                                    Y           = mplesetup$response, 
                                    X           = mplesetup$predictor, 
                                    weights     = mplesetup$weights, 
                                    prior.mean  = prior.mean, 
                                    prior.sigma = prior.sigma, 
                                    V           = S.prop, 
                                    thin        = thin, 
                                    mcmc        = main.iters, 
                                    burnin      = burn.in, 
                                    calibr.info = calibr.info, 
                                    seed        = seed, 
                                    logfun      = TRUE))
  
  message(" > Model evidence estimation")
  
  T0 <- as.mcmc(T0[(burn.in + 1):main.iters, ])
  AR <- round(1 - rejectionRate(T0)[1], 2)
  
  names(AR)  <- NULL
  ess        <- round(effectiveSize(T0), 0)
  names(ess) <- model$coef.names
  theta.sum  <- summary(T0)
  thetahat   <- theta.sum$statistics[, "Mean"]
  
  log.post <- adjusted_logPP(thetahat, 
                             Y           = mplesetup$response, 
                             X           = mplesetup$predictor, 
                             weights     = mplesetup$weights, 
                             prior.mean  = prior.mean, 
                             prior.sigma = prior.sigma, 
                             calibr.info = calibr.info)
  
  g       <- sample(1:nrow(T0), num.samples, replace = TRUE)
  theta.g <- T0[g, ]
  
  q.g <- dmvnorm(theta.g, mean = thetahat, sigma = S.prop, log = FALSE)
  
  lik.star <- adjusted_logPL(thetahat, 
                             Y           = mplesetup$response, 
                             X           = mplesetup$predictor, 
                             weights     = mplesetup$weights, 
                             calibr.info = calibr.info)
  
  lik.g <- apply(theta.g, 1, function(i) {
    adjusted_logPL(i, 
                   Y           = mplesetup$response, 
                   X           = mplesetup$predictor, 
                   weights     = mplesetup$weights, 
                   calibr.info = calibr.info)
  })
  
  alpha.g <- sapply(lik.g, function(l) min(1, exp(lik.star - l)))
  theta.j <- rmvnorm(num.samples, mean = thetahat, sigma = S.prop)
  
  lik.j <- apply(theta.j, 1, function(i) {
    adjusted_logPL(i, Y = mplesetup$response, 
                   X           = mplesetup$predictor, 
                   weights     = mplesetup$weights, 
                   calibr.info = calibr.info)
  })
  
  alpha.j     <- sapply(lik.j, function(l) min(1, exp(l - lik.star)))
  pi.hat      <- mean(alpha.g * q.g)/mean(alpha.j)
  
  logEvidence <- calibr.info$logC + log.post - log(pi.hat)
  
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)
  
  out <- list(formula      = formula, 
              Theta        = T0, 
              AR           = AR, 
              ess          = ess, 
              log.evidence = logEvidence, 
              dim          = dim, 
              specs        = model$coef.names, 
              Time         = runtime)
  class(out) <- "bergm"
  return(out)
}