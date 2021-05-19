#' Calibrating misspecified Bayesian ERGMs
#'
#' Function to transform a sample from the pseudo-posterior 
#' to one that is approximately sampled from the intractable 
#' posterior distribution.
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
#' @param burn.in count; number of burn-in iterations at the beginning of an MCMC run for the pseudo-posterior estimation.
#'
#' @param main.iters count; number of MCMC iterations after burn-in for the pseudo-posterior estimation.
#'
#' @param aux.iters count; number of auxiliary iterations used for drawing the first network from the ERGM likelihood (Robbins-Monro). See \code{\link[ergm]{control.simulate.formula}}.
#'
#' @param V.proposal count; diagonal entry for the multivariate Normal proposal.
#' By default set to 1.5.
#'
#' @param thin count; thinning interval used in the simulation for the pseudo-posterior estimation. The number of MCMC iterations must be divisible by this value.
#'
#' @param rm.iters count; number of iterations for the Robbins-Monro stochastic approximation algorithm.
#' 
#' @param rm.a scalar; constant for sequence alpha_n (Robbins-Monro).
#' 
#' @param rm.alpha scalar; noise added to gradient (Robbins-Monro).
#' 
#' @param n.aux.draws count; number of auxiliary networks drawn from the ERGM likelihood (Robbins-Monro). See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param aux.thin count; number of auxiliary iterations between network draws after the first network is drawn (Robbins-Monro). See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param estimate If "MLE" (the default), then an approximate maximum likelihood estimator is used as a starting point in the Robbins-Monro algorithm. If "CD" , the Monte-Carlo contrastive divergence estimate is returned. See \code{\link[ergm]{ergm}}.
#' 
#' @param seed integer; seed for the random number generator. See \code{set.seed}.
#' 
#' @param ... Additional arguments, to be passed to the ergm function. See \code{\link[ergm]{ergm}}.
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
#'
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
                   estimate    = c("MLE","CD"),
                   seed        = 1,
                   ...){
  
  y     <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  sy    <- summary(formula)    
  dim   <- length(sy)
  
  if (is.null(prior.mean))  prior.mean  <- rep(0, dim)
  if (is.null(prior.sigma)) prior.sigma <- diag(100, dim, dim)
  
  control <- control.ergm(MCMC.burnin     = aux.iters,
                          MCMC.interval   = aux.thin,
                          MCMC.samplesize = n.aux.draws,
                          seed            = seed)
  
  proposal <- ergm_proposal(object      = ~., 
                            constraints = ~., 
                            arguments   = control$MCMC.prop.args, 
                            nw          = y)
  
  estimate <- match.arg(estimate)
  
  logpp_short <- function(theta,
                          Y,
                          X,
                          weights,
                          prior.mean,   
                          prior.sigma,
                          optimPL){
    xtheta    <- c(X %*% theta)
    log.like  <- sum(dbinom(weights * Y, weights, exp(xtheta) / (1 + exp(xtheta)), log = TRUE))
    log.prior <- dmvnorm(theta, mean = prior.mean, prior.sigma, log = TRUE) 
    if (optimPL == FALSE) out <- log.like + log.prior else out <- -(log.like + log.prior)
    return(out)
  }
  
  score_logPP <- function(theta,          
                          Y,
                          X,
                          weights,
                          prior.mean,
                          prior.sigma,
                          optimPL){
    
    score.log.prior <- -solve(prior.sigma, (theta - prior.mean))
    p               <- c( exp(as.matrix(X) %*% theta ) / ( 1 + exp(as.matrix(X) %*% theta) ) )
    deriv.logpl     <- c( t(X) %*% (weights*(Y - p)) )
    if (optimPL == TRUE) out <- deriv.logpl + score.log.prior
    return(out)
  }
  
  rm_ergm <- function(formula,
                      prior.mean,      
                      prior.sigma,
                      y,
                      model,
                      sy,
                      dim,
                      control,
                      proposal,
                      rm.iters,            
                      rm.a,                 
                      rm.alpha,          
                      init,            
                      aux.iters, 
                      n.aux.draws, 
                      aux.thin){     
    
    theta <- array(0, c(dim = rm.iters, dim))  
    theta[1, ] <- init
    
    for (i in 2:rm.iters) {
      
      # Addition:
      y0 <- simulate(formula, 
                     coef        = theta[i - 1,], 
                     nsim        = 1, 
                     control     = control.simulate(MCMC.burnin = 1,
                                                    MCMC.interval = 1),
                     return.args = "ergm_state")$object
      
      # Addition:
      z <- as.matrix( ergm_MCMC_sample(y0,
                            theta   = theta[i - 1,],
                            stats0  = sy,
                            control = control)$stats[[1]] )

      estgrad    <- -apply(z, 2, mean) - solve(prior.sigma, (theta[i - 1,] - prior.mean))
      theta[i, ] <- theta[i - 1, ]  + ((rm.a/i) * (rm.alpha + estgrad))
    }
    out <- list(Theta = theta) 
    return(out)
  }
  
  mplesetup <- ergmMPLE(formula)
  data.glm.initial <- cbind(mplesetup$response, 
                            mplesetup$weights, 
                            mplesetup$predictor)
  colnames(data.glm.initial) <- c("responses", 
                                  "weights", 
                                  colnames(mplesetup$predictor))
  
  Vcov.MPLE <- vcov(glm(mplesetup$response ~. - 1, 
                        data    = data.frame(mplesetup$predictor), 
                        weights = mplesetup$weights, 
                        family  = "binomial"))
  
  Sigma.proposal <- diag(rep(V.proposal, dim), dim, dim)
  S.prop <- Sigma.proposal %*% solve((solve(prior.sigma) + solve(Vcov.MPLE))) %*% Sigma.proposal 
  
  suppressMessages( mple <- ergm(formula, estimate = "MPLE", verbose = FALSE) )
  
  message(" > MCMC start")
  clock.start <- Sys.time() 
  
  capture.output(unadj.sample <- MCMCmetrop1R(logpp_short,
                                              theta.init  = mple$coefficients,
                                              Y           = mplesetup$response,
                                              X           = mplesetup$predictor,
                                              weights     = mplesetup$weights,
                                              optimPL     = FALSE,
                                              prior.mean  = prior.mean,
                                              prior.sigma = prior.sigma,
                                              V           = S.prop,
                                              thin        = thin, 
                                              mcmc        = main.iters, 
                                              burnin      = burn.in,
                                              seed        = seed,
                                              logfun      = TRUE
                                              ))
  
  message(" > MAP estimation")
  
  suppressMessages(rob.mon.init <- ergm(formula, estimate = estimate, verbose = FALSE, control = control.ergm(seed = seed), ...))
  
  theta.star <- rm_ergm(formula, 
                        rm.iters    = rm.iters,
                        rm.a        = rm.a,
                        rm.alpha    = rm.alpha,
                        init        = rob.mon.init$coefficients,
                        aux.iters   = aux.iters,
                        n.aux.draws = n.aux.draws, 
                        aux.thin    = aux.thin,  
                        prior.mean  = prior.mean,     
                        prior.sigma = prior.sigma,
                        y           = y,
                        model       = model,
                        sy          = sy,
                        dim         = dim,
                        control     = control,
                        proposal    = proposal
  )
  
  theta.PLstar <- optim(par         = summary(unadj.sample)$statistics[,1],
                        fn          = logpp_short,
                        gr          = score_logPP,
                        lower       = rob.mon.init$coefficients - 3,
                        upper       = rob.mon.init$coefficients + 3,
                        Y           = mplesetup$response,
                        X           = mplesetup$predictor,
                        weights     = mplesetup$weights,
                        optimPL     = TRUE,
                        prior.mean  = prior.mean,     
                        prior.sigma = prior.sigma,
                        hessian     = TRUE,
                        method      = "L-BFGS-B")
  
  message(" > Curvature Adjustment")

  y0 <- simulate(formula, 
                 coef        = theta.star$Theta[rm.iters, ], 
                 nsim        = 1, 
                 control     = control.simulate(MCMC.burnin = 1,
                                                MCMC.interval = 1),
                 return.args = "ergm_state")$object
  
  z <- as.matrix( ergm_MCMC_sample(y0,
                                   theta   = theta.star$Theta[rm.iters, ],
                                   stats0  = sy,
                                   control = control)$stats[[1]] )
  sim.samples <- sweep(z, 2, sy, '+')
  
  Hessian.post.truelike <- -cov(sim.samples) - solve(prior.sigma)  
  
  chol.true.posterior <- chol(-Hessian.post.truelike)
  chol.PL.posterior   <- chol(-theta.PLstar$hessian) 
  
  W <- solve(chol.PL.posterior) %*% chol.true.posterior
  V <- solve(W)
  
  corrected.sample <- t(apply(data.frame(unadj.sample), 1, 
                              function(x) {
                                c(V %*% (unlist(x) - theta.PLstar$par ) + 
                                    theta.star$Theta[rm.iters, ])}))
  
  clock.end <- Sys.time()
  runtime   <- difftime(clock.end, clock.start)  
  
  FF  <- mcmc(corrected.sample)
  ess <- round(effectiveSize(FF), 0)
  names(ess) <- model$coef.names
  
  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL
  
  out <- list(Theta_star = theta.star$Theta[rm.iters, ],
              Theta      = FF,
              Time       = runtime,
              formula    = formula,
              AR         = AR,
              ess        = ess,
              dim        = dim,
              specs      = model$coef.names)
  
  class(out) <- "bergm"
  return(out)
}