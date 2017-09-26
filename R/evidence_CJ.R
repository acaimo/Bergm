#' Marginal likelihood estimation 
#'
#' Function to estimate the marginal likelihood with Chib and Jeliazkov's method, 
#' based on the adjusted pseudolikelihood function.
#' 
#' @param ergm.formula formula; an \code{R} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link{network}} object
#' and <model terms> are \code{\link{ergm-terms}}.
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
#' Bouranis, L., Friel, N., and Maire, F. (2017). Bayesian model selection for exponential random graph models via
#' adjusted pseudolikelihoods. \url{https://arxiv.org/abs/1706.06344}
#'
#' @examples
#' # Load the florentine marriage network:
#' data(florentine)
#'
#' flo.formula <- flomarriage ~ edges + kstar(2)
#'
#' info.adjustPL <- adjustPL(ergm.formula = flo.formula,
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
#' Chib.est.evidence <- evidence_CJ(ergm.formula= flo.formula,
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
#'                                                                                                     
#' @import network
#' @import ergm
#' @import mvtnorm
#' @import MCMCpack
#' @import Matrix
#'
#' @export
#'
evidence_CJ <- function(ergm.formula,
                          prior.mean,   
                          prior.sigma,
                          nits,
                          burnin,
                          thin        = 1,
                          num.samples = 5000,   # integer; number of samples used in the estimate
                          tunePL      = 2,
                          seed        = NA,
                          calibr.info = NULL)
  { 
  
  #==========================
  # Sub-routines:
  mspecs <- function(x) {
      y       <- ergm.getnetwork(x)
      model   <- ergm.getmodel(x, y)
      model$coef.names
  }# End function
  
  expit <- function(x) exp(x)/(1+exp(x)) 
  
  logPL <- function(theta,
                    y,
                    X,
                    weights){
    xtheta  <- c(X%*%theta)
    
    p       <- expit(xtheta)
    log.like<- sum( dbinom(weights*y, weights, expit(xtheta), log = TRUE) ) 
    
    return(log.like)
  }# End function
  
  logPL.corr <- function(theta,
                         y,
                         X,
                         weights,
                         calibr.info){
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL )
    xtheta  <- c(X%*%theta_transf)
    p       <- expit(xtheta)
    log.like<- sum( dbinom(weights*y, weights, expit(xtheta), log = TRUE) )
    
    return(log.like)
  }# End function
  
  logPP <- function(theta,
                    y,
                    X,
                    weights,
                    prior.mean,   
                    prior.sigma){
    
    xtheta   <- c(X%*%theta)
    p        <- expit(xtheta)
    log.like <- sum( dbinom(weights*y, weights, expit(xtheta), log = TRUE))
    log.prior<- dmvnorm(theta, mean = prior.mean, prior.sigma, log = TRUE) 
    
    return(log.like + log.prior)
  }# End function
  
  logPP.corr <- function(theta,
                         y,
                         X,
                         weights,
                         prior.mean,   
                         prior.sigma,
                         calibr.info){
    
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL )
    
    xtheta   <- c(X%*%theta_transf)
    p        <- expit(xtheta)
    log.like <- sum( dbinom(weights*y, weights, expit(xtheta), log = TRUE) )
    
    log.prior<- dmvnorm(theta, mean = prior.mean, prior.sigma, log = TRUE) 
    
    out <- log.like + log.prior
    return(out)
    
  }# End function
  
  #==========================
  Hessian.corrected.logpl <- function(theta,          
                                      y,
                                      X,
                                      weights,
                                      calibr.info){
    
    p <- exp(as.matrix(X)%*%theta)/ ( 1+exp(as.matrix(X)%*%theta) )
    W <- Diagonal( x=as.vector( weights*p*(1-p) )) 
    
    Hessian <- - t( as.matrix(X) ) %*% W %*% as.matrix(X)
    Hessian <- as.matrix(Hessian)
    
    return( t(calibr.info$W) %*%Hessian%*% calibr.info$W ) 
  }# End function
  
  #==========================
  mplesetup                 <- ergmMPLE(ergm.formula)
  data.glm.initial          <- cbind(mplesetup$response, 
                                     mplesetup$weights, 
                                     mplesetup$predictor)
  colnames(data.glm.initial)<- c("responses", 
                                 "weights", 
                                 colnames(mplesetup$predictor))
  
  Vcov.MPLE <- vcov(glm(mplesetup$response ~ . - 1, 
                        data    = data.frame(mplesetup$predictor), 
                        weights = mplesetup$weights, 
                        family  ="binomial"))
  
  mdims <- length(mspecs(ergm.formula))
  
  tune.mat.PL      <- diag(rep(tunePL, mdims))
  B0               <- solve(prior.sigma)      
  
  #==========================
  clock.start <- Sys.time() 
  
  #==========================
  H <- Hessian.corrected.logpl(calibr.info$Theta_PL,          
                               y          = data.glm.initial[,"responses"],
                               X          = data.glm.initial[,-c(1,2)],
                               weights    = data.glm.initial[,"weights"],
                               calibr.info= calibr.info)
  
  PL.prop.sigma.mat<- tune.mat.PL%*%solve( B0 - H )%*%tune.mat.PL 
  
  #==========================
  message("---MCMC start---")
  
  T0 <- MCMCmetrop1R(logPP.corr,
                     theta.init = calibr.info$Theta_PL,
                     y          = mplesetup$response,
                     X          = mplesetup$predictor,
                     weights    = mplesetup$weights,
                     prior.mean = prior.mean,
                     prior.sigma= prior.sigma,
                     V          = PL.prop.sigma.mat,
                     thin       = thin, 
                     mcmc       = nits,
                     burnin     = burnin, 
                     calibr.info= calibr.info,
                     seed       = seed,
                     logfun     = TRUE)
  
  message("---MCMC end---")
  
  #==========================
  # Info about the chain:
  T0          <- as.mcmc( T0[(burnin+1):nits,]) 
  AR          <- round(1-rejectionRate(T0)[1], 2)
  names(AR)   <- NULL
  ess         <- round(effectiveSize(T0),0)
  names(ess)  <- mspecs(ergm.formula)

  theta.sum  <- summary(T0)
  thetahat   <- theta.sum$statistics[,"Mean"]
  
  message("---Evidence estimation---")
  
  log.post <- logPP.corr(thetahat,
                         y          = mplesetup$response,
                         X          = mplesetup$predictor,
                         weights    = mplesetup$weights,
                         prior.mean = prior.mean,
                         prior.sigma= prior.sigma,
                         calibr.info= calibr.info) 
  
  g       <- sample(1:nrow(T0), num.samples, replace=TRUE)
  theta.g <- T0[g,]
  
  q.g <- dmvnorm(theta.g,
                 mean  = thetahat,
                 sigma = PL.prop.sigma.mat,
                 log   = FALSE) 
  
  lik.star<- logPL.corr(thetahat,
                        y          = mplesetup$response,
                        X          = mplesetup$predictor,
                        weights    = mplesetup$weights,
                        calibr.info= calibr.info)
  
  lik.g <- apply(theta.g, 1, 
                 function (i) {
                   logPL.corr(i,
                              y          = mplesetup$response,
                              X          = mplesetup$predictor,
                              weights    = mplesetup$weights,
                              calibr.info= calibr.info)
                 }) 
  
  alpha.g <- sapply(lik.g, function(l) min(1, exp(lik.star-l)) )
  
  theta.j <- rmvnorm(num.samples, mean=thetahat, sigma=PL.prop.sigma.mat)
  lik.j   <- apply(theta.j, 1, 
                   function (i) {
                     logPL.corr(i,
                                y          = mplesetup$response,
                                X          = mplesetup$predictor,
                                weights    = mplesetup$weights,
                                calibr.info= calibr.info) 
                   }) 
  
  alpha.j <- sapply(lik.j, function(l) min(1,exp(l-lik.star)))
  pi.hat  <- mean(alpha.g*q.g)/mean(alpha.j)
  
  logEvidence <- log(calibr.info$C) + log.post - log(pi.hat)
  
  #==========================
  clock.end <- Sys.time()
  runtime   <- difftime(clock.end, clock.start)  

  #==========================
  out <- list(Theta       = T0,
              dim         = mdims,
              log.evidence= logEvidence,
              formula     = ergm.formula,
              AR          = AR,
              ESS         = ess,
              dim         = mdims,
              specs       = mspecs(ergm.formula),
              Time        = runtime)
  class(out) <- "evidence_CJ"
  
  return(out)
  
}# End function
