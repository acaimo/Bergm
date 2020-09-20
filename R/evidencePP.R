#' Evidence estimation via power posteriors
#'
#' Function to estimate the evidence (marginal likelihood) with Power posteriors, 
#' based on the adjusted pseudolikelihood function.
#' 
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}.
#' 
#' @param prior.mean vector; mean vector of the multivariate Normal prior.
#' By default set to a vector of 0's.
#' 
#' @param prior.sigma square matrix; variance/covariance matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100.
#' 
#' @param aux.iters count; number of auxiliary iterations used for drawing the first network from the ERGM likelihood. See \code{\link[ergm]{control.simulate.formula}} and \code{\link[Bergm]{ergmAPL}}.
#' 
#' @param n.aux.draws count; number of auxiliary networks drawn from the ERGM likelihood. See \code{\link[ergm]{control.simulate.formula}} and \code{\link[Bergm]{ergmAPL}}.
#' 
#' @param aux.thin count; number of auxiliary iterations between network draws after the first network is drawn. See \code{\link[ergm]{control.simulate.formula}} and \code{\link[Bergm]{ergmAPL}}.
#' 
#' @param ladder count; length of temperature ladder (>=3). See \code{\link[Bergm]{ergmAPL}}.
#' 
#' @param main.iters count; number of MCMC iterations after burn-in for the adjusted pseudo-posterior estimation.
#' 
#' @param burn.in count; number of burn-in iterations at the beginning of an MCMC run for the adjusted pseudo-posterior estimation.
#' 
#' @param thin count; thinning interval used in the simulation for the adjusted pseudo-posterior estimation. The number of MCMC iterations must be divisible by this value.
#' 
#' @param V.proposal count; diagonal entry for the multivariate Normal proposal.
#' By default set to 1.5.
#' 
#' @param seed integer; seed for the random number generator. 
#' See \code{set.seed} and \code{\link[MCMCpack]{MCMCmetrop1R}}. 
#' 
#' @param temps numeric vector; inverse temperature ladder, \eqn{t \in [0,1]}.
#'
#' @param estimate If "MLE" (the default), then an approximate maximum likelihood estimator is returned. If "CD" , the Monte-Carlo contrastive divergence estimate is returned. See \code{\link[ergm]{ergm}}.
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
#' @examples
#' \dontrun{
#' # Load the florentine marriage network:
#' data(florentine)
#' 
#' PPE <- evidencePP(flomarriage ~ edges + kstar(2),
#'                   aux.iters   = 500, 
#'                   noisy.nsim  = 50,   
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
#'
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
                       ...) 
{
  y     <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  sy    <- summary(formula)
  dim   <- length(sy)
  
  if (dim == 1) stop("Model dimension must be greater than 1")
  if (is.null(prior.mean)) prior.mean <- rep(0, dim)
  if (is.null(prior.sigma)) prior.sigma <- diag(100, dim, dim)
  if (is.null(temps)) temps <- seq(0, 1, length.out = 50)^5
  
  estimate <- match.arg(estimate)
  
  expit <- function(x) exp(x)/(1 + exp(x))
  
  adjusted_logPL <- function(theta, 
                             Y, 
                             X, 
                             weights, 
                             calibr.info) {
    theta_transf <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL)
    xtheta   <- c(X %*% theta_transf)
    log.like <- sum(dbinom(weights * Y, weights, expit(xtheta), log = TRUE))
    return(log.like)
  }
  
  adjusted_logPP <- function(theta, 
                             Y, 
                             X, 
                             weights, 
                             prior.mean, 
                             prior.sigma, 
                             temperature, 
                             calibr.info) {
    theta_transf <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL)
    xtheta <- c(X %*% theta_transf)
    p <- expit(xtheta)
    log.like <- sum(dbinom(weights * Y, weights, expit(xtheta), log = TRUE))
    log.prior <- dmvnorm(theta, mean = prior.mean, prior.sigma, log = TRUE)
    return(temperature * log.like + log.prior)
  }
  
  score_logPPt <- function(theta, 
                           Y, 
                           X, 
                           weights, 
                           prior.mean, 
                           prior.sigma, 
                           calibr.info, 
                           temperature) {
    score.log.prior <- -solve(prior.sigma, (theta - prior.mean))
    theta_transf    <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL)
    p <- c(exp(as.matrix(X) %*% theta_transf)/(1 + exp(as.matrix(X) %*% theta_transf)))
    deriv.logpl <- c(t(X) %*% (weights * (Y - p)))
    out <- temperature * (t(calibr.info$W) %*% deriv.logpl) + score.log.prior
    return(out)
  }
  
  Hessian_adjusted_logPL <- function(theta, 
                                     X, 
                                     weights, 
                                     calibr.info) {
    p <- exp(as.matrix(X) %*% theta)/(1 + exp(as.matrix(X) %*% theta))
    W <- Diagonal(x = as.vector(weights * p * (1 - p)))
    Hessian <- -t(as.matrix(X)) %*% W %*% as.matrix(X)
    Hessian <- as.matrix(Hessian)
    return(t(calibr.info$W) %*% Hessian %*% calibr.info$W)
  }
  
  numtemp     <- length(temps)
  pplist      <- list()
  acceptances <- rep(0, numtemp)
  cv.ell.D2   <- rep(0, numtemp)
  cv.vll.D2   <- rep(0, numtemp)
  htheta.D2   <- list()
  
  clock.start <- Sys.time()
  
  message(" > Pseudolikelihood adjustment")
  info.adjustPL <- ergmAPL(formula     = formula, 
                           aux.iters   = aux.iters, 
                           n.aux.draws = n.aux.draws, 
                           aux.thin    = aux.thin, 
                           ladder      = ladder,
                           estimate    = estimate,
                           seed        = seed,
                           ...)
  
  calibr.info <- list(Theta_MLE = info.adjustPL$Theta_MLE, 
                      Theta_PL  = info.adjustPL$Theta_PL, 
                      W         = info.adjustPL$W, 
                      logC      = info.adjustPL$logC)
  
  mplesetup        <- ergmMPLE(formula)
  data.glm.initial <- cbind(mplesetup$response, mplesetup$weights, mplesetup$predictor)
  colnames(data.glm.initial) <- c("responses", "weights", colnames(mplesetup$predictor))
  Vcov.MPLE <- vcov(glm(mplesetup$response ~ . - 1, 
                        data    = data.frame(mplesetup$predictor), 
                        weights = mplesetup$weights, 
                        family  = "binomial"))
  
  Sigma.proposal <- diag(rep(V.proposal, dim), dim, dim)
  H <- Hessian_adjusted_logPL(theta       = calibr.info$Theta_PL, 
                              X           = data.glm.initial[, -c(1, 2)], 
                              weights     = data.glm.initial[,"weights"], 
                              calibr.info = calibr.info)
  
  S.prop <- Sigma.proposal %*% solve(solve(prior.sigma) - H) %*% Sigma.proposal
  rownames(S.prop) <- colnames(S.prop) <- rownames(Vcov.MPLE)
  taup        <- 1/sqrt(diag(S.prop))
  tau0        <- solve(prior.sigma)
  a.tilde     <- log(tau0[1, 1]/taup)/log(temps[2])
  covar.temps <- list()
  
  for (i in numtemp:1) {
    covar.temps[[i]] <- S.prop
    diag(covar.temps[[i]]) <- diag(S.prop)/(temps[i]^a.tilde[1])
  }
  covar.temps[[1]] <- covar.temps[[2]]
  
  message(" > MCMC run")
  
  for (l in numtemp:1) {
    
    if (l == 1) {
      
      #set.seed(seed)
      
      pplist[[l]]    <- rmvnorm(main.iters, mean = prior.mean, sigma = prior.sigma)
      acceptances[l] <- 1
      
    } else {
      if (l == numtemp) l0 <- calibr.info$Theta_PL else l0 <- c(tail(red)[6, ])
      
      capture.output(red <- MCMCmetrop1R(adjusted_logPP, 
                                         mcmc        = main.iters, 
                                         burnin      = burn.in, 
                                         thin        = thin, 
                                         theta.init  = l0, 
                                         seed        = NA, 
                                         V           = covar.temps[[l]], 
                                         Y           = mplesetup$response, 
                                         X           = mplesetup$predictor, 
                                         weights     = mplesetup$weights, 
                                         calibr.info = calibr.info, 
                                         prior.mean  = prior.mean, 
                                         prior.sigma = prior.sigma, 
                                         temperature = temps[l], 
                                         logfun      = TRUE))
      pplist[[l]] <- red
      acceptances[l] <- round(1 - rejectionRate(red)[1], 3)
    }
  }
  
  message(" > Model evidence estimation")
  
  names(acceptances) <- paste("Temp", "_", 1:numtemp, " = ", round(temps, 5), sep = "")
  grad.logpp.list <- logPL.list <- list()
  
  for (j in 1:numtemp) {
    if (j == 1) {
      grad.logpp.list[[j]] <- t(apply(data.frame(pplist[[j]]), 1, function(x) -solve(prior.sigma, (x - prior.mean))))
      logPL.list[[j]] <- apply(pplist[[j]], 1, function(x) dmvnorm(x,prior.mean, prior.sigma, log = TRUE))
    } else {
      grad.logpp.list[[j]] <- t(apply(data.frame(pplist[[j]]), 
                                      1, function(x) {
                                        score_logPPt(theta       = x, 
                                                     Y           = mplesetup$response, 
                                                     X           = mplesetup$predictor, 
                                                     weights     = mplesetup$weights, 
                                                     prior.mean  = prior.mean, 
                                                     prior.sigma = prior.sigma, 
                                                     temperature = temps[j], 
                                                     calibr.info = calibr.info)
                                      }))
      logPL.list[[j]] <- apply(data.frame(pplist[[j]]), 
                               1, function(x) {
                                 adjusted_logPL(theta       = x, 
                                                Y           = mplesetup$response, 
                                                X           = mplesetup$predictor, 
                                                weights     = mplesetup$weights, 
                                                calibr.info = calibr.info)
                               })
    }
    k <- dim * (dim + 3)/2
    l <- 2 * dim + 1
    w.mat <- matrix(0, nrow = main.iters, ncol = k)
    w.mat[, 1:dim] <- grad.logpp.list[[j]]
    w.mat[, (dim + 1):(2 * dim)] <- as.matrix(pplist[[j]]) * grad.logpp.list[[j]] + rep(1, main.iters)
    
    for (k1 in (1:(dim - 1))) {
      for (k2 in (k1 + 1):dim) {
        w.mat[, l] <- as.matrix(pplist[[j]])[, k1] * 
          grad.logpp.list[[j]][, k2] + 
          as.matrix(pplist[[j]])[,k2] * grad.logpp.list[[j]][, k1]
        l <- l + 1
      }
    }
    phi.D2         <- c(-solve(var(w.mat)) %*% c(cov(w.mat, logPL.list[[j]])))
    htheta.D2[[j]] <- c(phi.D2 %*% t(w.mat))
    cv.ell.D2[j]   <- mean(logPL.list[[j]] + htheta.D2[[j]])
    cv.vll.D2[j]   <- var(logPL.list[[j]] + htheta.D2[[j]])
  }
  
  clock.end <- Sys.time()
  runtime   <- difftime(clock.end, clock.start)
  
  ppml <- function(cv.ell.D2, cv.vll.D2, tempvec) {
    N <- length(cv.ell.D2)
    cv.res.D2.mts <- 0
    for (i in 1:(N - 1)) {
      wts <- tempvec[i + 1] - tempvec[i]
      cv.res.D2.mts <- cv.res.D2.mts + wts * ((cv.ell.D2[i +1] + cv.ell.D2[i])/2) - 
                       ((wts^2)/12) * (cv.vll.D2[i + 1] - cv.vll.D2[i])
    }
    return(cv.res.D2.mts)
  }
  pp.estimates <- ppml(cv.ell.D2 = cv.ell.D2, cv.vll.D2 = cv.vll.D2, tempvec = temps)
  
  log.evidence <- calibr.info$logC + pp.estimates
  ess          <- round(effectiveSize(pplist[[numtemp]]), 0)
  names(ess)   <- model$coef.names
  
  out <- list(formula = formula, 
              Theta   = pplist[[numtemp]], 
              AR_temp = acceptances, 
              AR      = acceptances[numtemp], 
              ess     = ess, 
              log.evidence = log.evidence, 
              dim     = dim, 
              specs   = model$coef.names, 
              Time    = runtime)
  class(out) <- "bergm"
  return(out)
}