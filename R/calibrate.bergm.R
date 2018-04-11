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
#' @param iters count; Iterations for the Robbins-Monro stochastic approximation algorithm.
#' 
#' @param a scalar; Constant for sequence alpha_n (Robbins-Monro).
#' 
#' @param alpha scalar; Noise added to gradient (Robbins-Monro).
#' 
#' @param aux.iters count; Number of proposals before any MCMC sampling is done (Robbins-Monro). See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param noisy.nsim count; Number of TNT draws (Robbins-Monro). See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param noisy.thin count; Number of proposals between sampled statistics (Robbins-Monro). See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param prior.mean vector; Prior means.
#' 
#' @param prior.sigma matrix; Prior covariance matrix.
#' 
#' @param thin count; Thinning interval used in the simulation for the pseudo-posterior estimation. The number of MCMC iterations must be divisible by this value.
#' 
#' @param mcmc count; Number of MCMC iterations after burn-in for the pseudo-posterior estimation.
#' 
#' @param burnin count; Number of burn-in iterations at the beginning of an MCMC run for the pseudo-posterior estimation.
#' 
#' @param tunePL count; Tuning parameter for the Metropolis sampling for the pseudo-posterior estimation.
#'
#' @references
#' Bouranis, L., Friel, N., and Maire, F. (2015). Bayesian inference for misspecified exponential
#' random graph models. \url{https://arxiv.org/abs/1510.00934}
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data(florentine)
#'                                  
#' # Calibrated pseudo-posterior:
#'
#' cpp.flo <- calibrate.bergm(flomarriage ~ edges + kstar(2),
#'                            aux.iters = 3000,
#'                            mcmc = 10000,  
#'                            burnin = 500,
#'                            tunePL = 2)
#'                                                    
#' # MCMC diagnostics and posterior summaries:
#' 
#' bergm.output(cpp.flo)
#' 
#' # Bayesian goodness-of-fit test:
#' 
#' bgof(cpp.flo,
#'      aux.iters = 500,
#'      sample.size = 50,
#'      n.deg = 10,
#'      n.dist = 9,
#'      n.esp = 6)
#'}
#'
#' @export
#'

calibrate.bergm<- function(formula,      
                           iters       = 500,            
                           a           = 0.001,                
                           alpha       = 0,          
                           aux.iters   = 5000, 
                           noisy.nsim  = 400,  
                           noisy.thin  = 50,   
                           prior.mean  = NULL,        
                           prior.sigma = NULL,      
                           thin        = 1,     
                           mcmc        = 4e04,  
                           burnin      = 1e04,  
                           tunePL      = 1){

  #=================================
  # SUB-ROUTINES
  
  # Obtain coefficient names for each model:
  mspecs <- function(x) {
    y <- ergm.getnetwork(x)
    model <- ergm_model(x, y)
    model$coef.names
  }
  
  # Log pseudo-posterior:
  expit <- function(x) exp(x) / (1 + exp(x)) 
  log_pseudo_post_short <- function(theta,
                                    y,
                                    X,
                                    weights,
                                    prior.mean,   
                                    prior.sigma,
                                    optimPL){
    xtheta <- c(X %*% theta)
    #p <- expit(xtheta) # useless?
    log.like <- sum(dbinom(weights * y, weights, expit(xtheta), log = TRUE))
    log.prior <- dmvnorm(theta, mean = prior.mean, prior.sigma, log = TRUE) 
    if (optimPL == FALSE) out <- log.like + log.prior else out <- -(log.like + log.prior)
    return(out)
  }
  
  # Score of the log pseudo-posterior:
  score.logpp <- function(theta,          
                          y,
                          X,
                          weights,
                          prior.mean,
                          prior.sigma,
                          optimPL){
    
    score.log.prior <- -solve(prior.sigma, (theta - prior.mean))
    
    p <- c( exp(as.matrix(X) %*% theta ) / ( 1 + exp(as.matrix(X) %*% theta) ) )
    deriv.logpl <- c( t(X) %*% (weights*(y - p))  )
    
    if (optimPL == TRUE) out <- deriv.logpl + score.log.prior
    return(out)
  }
  
  # Robins-Monro stochastic algorithm for finding
  # the MAP of the true posterior.
  robinsmon.ergm <- function(formula,           
                             iters,            
                             a,                 
                             alpha,          
                             init,            
                             aux.iters, 
                             noisy.nsim, 
                             noisy.thin,  
                             prior.mean,      
                             prior.sigma){     
    
    # --- Parameters for ergm.mcmcslave():
    # y       <- ergm.getnetwork(formula)
    # model   <- ergm.getmodel(formula, y)
    # Clist   <- ergm.Cprepare(y, model)
    # control <- control.simulate.formula(MCMC.burnin   = aux.iters, 
    #                                    MCMC.interval = noisy.thin)
    #
    # control$MCMC.samplesize <- noisy.nsim
    #
    # MHproposal <- MHproposal.ergm(object      = model, 
    #                              constraints = ~., 
    #                              arguments   = control$MCMC.prop.args, 
    #                              nw          = y, 
    #                              weights     = control$MCMC.prop.weights, 
    #                              class       = "c",
    #                              reference   = ~Bernoulli, 
    #                              response    = NULL)
    # ---

    sy <- summary(formula)    
    dim <- length(sy)
    theta <- array(0, c(dim = iters, dim))  
    theta[1, ] <- init
    
    pb <- txtProgressBar(min = 0, max = iters, style = 3)
    
    for(i in 2:iters){
      
      # ---
      #z <- ergm.mcmcslave(Clist      = Clist, 
      #                    MHproposal = MHproposal, 
      #                    eta0       = theta[i - 1,], 
      #                    control    = control, 
      #                    verbose    = FALSE)$s
      # ---
      sy1 <- simulate(formula, 
                      coef = theta[i - 1,],
                      nsim = noisy.nsim,
                      statsonly = TRUE,
                      control = control.simulate(MCMC.burnin = aux.iters, 
                                                 MCMC.interval = noisy.thin))
      
      z <- sweep(sy1, 2, sy, `-`)
      estgrad   <- -apply(z, 2, mean) - solve(prior.sigma, (theta[i-1,] - prior.mean))
      theta[i, ] <- theta[i - 1, ]  + ((a/i)*(alpha + estgrad))
      
      setTxtProgressBar(pb, i)
    }
    
    out <- list(Theta = theta) 
    return(out)
  }
  
  # Default prior specification:
  if (is.null(prior.mean)) prior.mean  <- rep(0, dim)
  if (is.null(prior.sigma)) prior.sigma <- diag(100, dim)
  
  # Get data in aggregated format:
  mplesetup <- ergmMPLE(formula)
  data.glm.initial <- cbind(mplesetup$response, 
                            mplesetup$weights, 
                            mplesetup$predictor)
  colnames(data.glm.initial) <- c("responses", 
                                  "weights", 
                                  colnames(mplesetup$predictor))
  
  # Variance-covariance matrix from MPLE:
  Vcov.MPLE <- vcov(glm(mplesetup$response ~. - 1, 
                        data    = data.frame(mplesetup$predictor), 
                        weights = mplesetup$weights, 
                        family  = "binomial"))

  
  # Tuning the variance-covariance matrix of the proposal distribution:
  tune.mat.PL <- diag(rep(tunePL, dim))
  B0 <- solve(prior.sigma)      
  PL.prop.sigma.mat<- tune.mat.PL %*% solve((B0 + solve(Vcov.MPLE))) %*% tune.mat.PL 
  
  # Obtain the MPLE for the model of interest:
  capture.output(mple <- ergm(formula, estimate = "MPLE") )
  
  message("---MCMC start---")
  clock.start <- Sys.time() 
  
  #==========================
  capture.output(unadj.sample <- MCMCmetrop1R(log_pseudo_post_short,
                                              theta.init  = mple$coef,
                                              y           = mplesetup$response,
                                              X           = mplesetup$predictor,
                                              weights     = mplesetup$weights,
                                              optimPL     = F,
                                              prior.mean  = prior.mean,
                                              prior.sigma = prior.sigma,
                                              V           = PL.prop.sigma.mat,
                                              thin        = thin, 
                                              mcmc        = mcmc, 
                                              burnin      = burnin,
                                              logfun      = TRUE))
  
  message("---MAP estimation---")
  
  #==========================
  # Use the MPLE to initialise the Robbins-Monro algorithm:
  capture.output(rob.mon.init <- ergm(formula, 
                                      estimate = "MLE"),
                                      control = control.ergm(MCMC.samplesize = 1000)) 
  sy <- summary(formula)
  
  #==========================
  # Estimate the MAP of the true posterior:
  theta.star <- robinsmon.ergm(formula, 
                               iters       = iters,
                               a           = a,
                               alpha       = alpha,
                               init        = rob.mon.init$coef,
                               aux.iters   = aux.iters,
                               noisy.nsim  = noisy.nsim, 
                               noisy.thin  = noisy.thin,  
                               prior.mean  = prior.mean,     
                               prior.sigma = prior.sigma)
  
  theta.PLstar <- optim(par         = summary(unadj.sample)$statistics[,1],
                        fn          = log_pseudo_post_short,
                        gr          = score.logpp,
                        lower       = rob.mon.init$coef-3,
                        upper       = rob.mon.init$coef+3,
                        y           = mplesetup$response,
                        X           = mplesetup$predictor,
                        weights     = mplesetup$weights,
                        optimPL     = TRUE,
                        prior.mean  = prior.mean,     
                        prior.sigma = prior.sigma,
                        hessian     = TRUE,
                        method      = "L-BFGS-B")
  
  message("---Curvature Adjustment---")
  
  #==========================
  # Simulate graphs from \theta^{*}::
  sim.samples <- simulate(formula, 
                          coef = theta.star$Theta[iters, ],
                          statsonly = TRUE,
                          nsim = noisy.nsim,
                          control = control.simulate(MCMC.burnin = aux.iters, 
                                                     MCMC.interval = noisy.thin))

  #==========================
  # Hessian of true log-posterior: 
  Hessian.post.truelike <- -cov(sim.samples) - solve(prior.sigma)  
  
  #==========================
  chol.true.posterior <- chol(-Hessian.post.truelike)
  chol.PL.posterior <- chol(-theta.PLstar$hessian) 
  
  #==========================
  # Calculate transformation matrix W:
  W <- solve(chol.PL.posterior) %*% chol.true.posterior
  V <- solve(W)
  
  #==========================
  # Correct the sample:
  corrected.sample <- t(apply(data.frame(unadj.sample), 1, 
                              function(x) {
                                c( V%*% (unlist(x) - theta.PLstar$par ) + 
                                          theta.star$Theta[iters, ] )
                              }
                             ) 
                       )
  clock.end <- Sys.time()
  runtime   <- difftime(clock.end, clock.start)  
  
  #==========================
  # Summarize the MCMC sample:
  corrected.sample.mcmc <- mcmc(corrected.sample)
  ess <- round(effectiveSize(corrected.sample.mcmc),0)
  names(ess) <- mspecs(formula)
  
  AR <- round(1 - rejectionRate(corrected.sample.mcmc)[1], 2)
  names(AR) <- NULL
  
  out <- list(Theta_star = theta.star$Theta[iters, ],
              Theta_PL = theta.PLstar$par,        
              W = W,                       
              Theta = corrected.sample.mcmc,
              Time = runtime,
              formula = formula,
              AR = AR,
              ess = ess,
              dim = dim,
              specs = names(ess))
  
  class(out) <- "calibrate.bergm"
  return(out)
}








