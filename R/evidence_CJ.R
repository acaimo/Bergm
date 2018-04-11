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
#' @param prior.mean vector; Prior means.
#' 
#' @param prior.sigma matrix; Prior covariance matrix.
#' 
#' @param nits count; Number of MCMC iterations after burn-in for the adjusted pseudo-posterior estimation.
#' 
#' @param burnin count; Number of burn-in iterations at the beginning of an MCMC run for the adjusted pseudo-posterior estimation.
#' 
#' @param thin count; Thinning interval used in the simulation for the adjusted pseudo-posterior estimation. The number of MCMC iterations must be divisible by this value.
#' 
#' @param num.samples integer; number of samples used in the marginal likelihood estimate. Must be <=(nits-burnin).
#' 
#' @param tunePL count; Tuning parameter for the Metropolis sampling for the pseudo-posterior estimation.
#' 
#' @param seed The seed for the random number generator. See \code{\link[MCMCpack]{MCMCmetrop1R}}. 
#' 
#' @param calibr.info list; Transformation parameters for
#' adjusting the pseudolikelihood function. See \code{\link[Bergm]{adjustPL}}. 
#'
#' @references
#' Caimo, A., & Friel, N. (2013). Bayesian model selection for exponential random graph models. 
#' Social Networks, 35(1), 11-24. \url{https://arxiv.org/abs/1201.2337}
#' 
#' Bouranis, L., Friel, N., and Maire, F. (2017). Bayesian model selection for exponential random graph models via
#' adjusted pseudolikelihoods. \url{https://arxiv.org/abs/1706.06344}
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network:
#' data(florentine)
#'
#' flo.formula <- flomarriage ~ edges + kstar(2)
#'
#' info.adjustPL <- adjustPL(formula = flo.formula,
#'                           aux.iters    = 100, 
#'                           noisy.nsim   = 50,   
#'                           noisy.thin   = 50,   
#'                           ladder       = 30,   
#'                           estimate     = "MLE",
#'                           control      = control.ergm(MCMC.samplesize=2000))
#'
#' # Add the output into a list:
#' calibration.list <- list(Theta_MLE= info.adjustPL$Theta_MLE,
#'                          Theta_PL = info.adjustPL$Theta_PL, 
#'                          W        = info.adjustPL$W, 
#'                          C        = info.adjustPL$C)
#'                          
#' # Specify location and shape of prior distribution: 
#' mean.priors  <- rep(0,2)
#' sigma        <- 5
#' sigma.priors <- diag(sigma,2)          
#'                                                 
#' # MCMC sampling and evidence estimation:
#' Chib.est.evidence <- evidence_CJ(formula= flo.formula,
#'                                  prior.mean   = mean.priors,   
#'                                  prior.sigma  = sigma.priors,
#'                                  nits         = 30000,
#'                                  burnin       = 5000,
#'                                  thin         = 1,
#'                                  num.samples  = 25000,
#'                                  tunePL       = 2,
#'                                  calibr.info  = calibration.list,
#'                                  seed         = 1)                         
#'                                    
#' # MCMC diagnostics and posterior summaries:
#' bergm.output(Chib.est.evidence)
#'   
#' # Log-marginal likelihood estimate:             
#' flo.model.logevidence <- Chib.est.evidence$log.evidence
#'}
#'
#' @export
#'
evidence_CJ <- function(formula,
                        prior.mean,   
                        prior.sigma,
                        nits,
                        burnin,
                        thin        = 1,
                        num.samples = 5000,   # integer; number of samples used in the estimate
                        tunePL      = 2,
                        seed        = NA,
                        calibr.info = NULL){ 
  #==========================
  # SUB-ROUTINES

  mspecs <- function(x) {
      y <- ergm.getnetwork(x)
      model <- ergm_model(x, y)
      model$coef.names
  }
  
  expit <- function(x) exp(x) / (1 + exp(x)) 
  
  logPL <- function(theta,
                    y,
                    X,
                    weights){
    xtheta  <- c(X %*% theta)
    p <- expit(xtheta)
    log.like <- sum(dbinom(weights * y, weights, expit(xtheta), log = TRUE)) 
    return(log.like)
  }
  
  logPL.corr <- function(theta,
                         y,
                         X,
                         weights,
                         calibr.info){
    
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL )
    xtheta  <- c(X %*% theta_transf)
    #p <- expit(xtheta) # useless?
    log.like <- sum( dbinom(weights * y, weights, expit(xtheta), log = TRUE) )
    return(log.like)
  }
  
  logPP <- function(theta,
                    y,
                    X,
                    weights,
                    prior.mean,   
                    prior.sigma){
    xtheta   <- c(X %*% theta)
    #p        <- expit(xtheta) # useless?
    log.like <- sum( dbinom(weights*y, weights, expit(xtheta), log = TRUE))
    log.prior<- dmvnorm(theta, prior.mean, prior.sigma, log = TRUE) 
    return(log.like + log.prior)
  }
  
  logPP.corr <- function(theta,
                         y,
                         X,
                         weights,
                         prior.mean,   
                         prior.sigma,
                         calibr.info){
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL )
    xtheta <- c(X %*% theta_transf)
    #p <- expit(xtheta) #useless?
    log.like <- sum(dbinom(weights * y, weights, expit(xtheta), log = TRUE))
    log.prior <- dmvnorm(theta, prior.mean, prior.sigma, log = TRUE) 
    return(log.like + log.prior)
  }
  
  #==========================
  Hessian.corrected.logpl <- function(theta,          
                                      y,
                                      X,
                                      weights,
                                      calibr.info){
    
    p <- exp(as.matrix(X) %*% theta) / ( 1 + exp(as.matrix(X) %*% theta))
    W <- Diagonal(x = as.vector( weights * p * (1 - p))) 
    Hessian <- as.matrix(-t( as.matrix(X) ) %*% W %*% as.matrix(X))
    return(t(calibr.info$W) %*% Hessian %*% calibr.info$W) 
  }
  
  #==========================
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
                        family  ="binomial"))
  
  dim <- length(mspecs(formula))
  tune.mat.PL <- diag(rep(tunePL, dim))
  B0 <- solve(prior.sigma)      
  
  #==========================
  H <- Hessian.corrected.logpl(calibr.info$Theta_PL,          
                               y           = data.glm.initial[,"responses"],
                               X           = data.glm.initial[, -c(1, 2)],
                               weights     = data.glm.initial[,"weights"],
                               calibr.info = calibr.info)
  
  PL.prop.sigma.mat<- tune.mat.PL %*% solve(B0 - H) %*% tune.mat.PL 
  
  #==========================
  # Start clock:
  clock.start <- Sys.time() 
  message("---MCMC start---")
  
  T0 <- MCMCmetrop1R(logPP.corr,
                     theta.init  = calibr.info$Theta_PL,
                     y           = mplesetup$response,
                     X           = mplesetup$predictor,
                     weights     = mplesetup$weights,
                     prior.mean  = prior.mean,
                     prior.sigma = prior.sigma,
                     V           = PL.prop.sigma.mat,
                     thin        = thin, 
                     mcmc        = nits,
                     burnin      = burnin, 
                     calibr.info = calibr.info,
                     seed        = seed,
                     logfun      = TRUE)
  
  #==========================
  # Info about the chain:
  T0          <- as.mcmc( T0[(burnin + 1):nits, ]) 
  AR          <- round(1 - rejectionRate(T0)[1], 2)
  names(AR)   <- NULL
  ess         <- round(effectiveSize(T0), 0)
  names(ess)  <- mspecs(formula)
  theta.sum   <- summary(T0)
  thetahat    <- theta.sum$statistics[, "Mean"]
  #message("---Evidence estimation---")
  
  log.post <- logPP.corr(thetahat,
                         y           = mplesetup$response,
                         X           = mplesetup$predictor,
                         weights     = mplesetup$weights,
                         prior.mean  = prior.mean,
                         prior.sigma = prior.sigma,
                         calibr.info = calibr.info) 
  
  g <- sample(1:nrow(T0), num.samples, replace = TRUE)
  theta.g <- T0[g, ]
  
  q.g <- dmvnorm(theta.g,
                 mean  = thetahat,
                 sigma = PL.prop.sigma.mat,
                 log   = FALSE) 
  
  lik.star<- logPL.corr(thetahat,
                        y           = mplesetup$response,
                        X           = mplesetup$predictor,
                        weights     = mplesetup$weights,
                        calibr.info = calibr.info)
  
  lik.g <- apply(theta.g, 1, 
                 function (i) {
                   logPL.corr(i,
                              y           = mplesetup$response,
                              X           = mplesetup$predictor,
                              weights     = mplesetup$weights,
                              calibr.info = calibr.info)
                   }
                 ) 
  
  alpha.g <- sapply(lik.g, function(l) min(1, exp(lik.star-l)) )
  
  theta.j <- rmvnorm(num.samples, mean = thetahat, sigma = PL.prop.sigma.mat)
  lik.j   <- apply(theta.j, 1, 
                   function (i) {
                     logPL.corr(i,
                                y           = mplesetup$response,
                                X           = mplesetup$predictor,
                                weights     = mplesetup$weights,
                                calibr.info = calibr.info) 
                   }) 
  
  alpha.j <- sapply(lik.j, function(l) min(1, exp(l - lik.star)))
  pi.hat  <- mean(alpha.g * q.g) / mean(alpha.j)
  
  logEvidence <- log(calibr.info$C) + log.post - log(pi.hat)
  
  #==========================
  # Stop clock:
  clock.end <- Sys.time()
  runtime   <- difftime(clock.end, clock.start)  

  #==========================
  out <- list(Theta        = T0,
              dim          = dim,
              log.evidence = logEvidence,
              formula      = formula,
              AR           = AR,
              ESS          = ess,
              specs        = names(ess),
              Time         = runtime)
  class(out) <- "evidence_CJ"
  return(out)
}
