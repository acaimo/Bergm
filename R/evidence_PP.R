#' Evidence estimation via power posteriors
#'
#' Function to estimate the evidence (marginal likelihood) with Power posteriors, 
#' based on the adjusted pseudolikelihood function.
#' 
#' @param formula formula; an \code{R} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link{network}} object
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
#' @param tunePL count; Tuning parameter for the Metropolis sampling for the pseudo-posterior estimation. 
#' 
#' @param seed The seed for the random number generator. See \code{\link[MCMCpack]{MCMCmetrop1R}}. 
#' 
#' @param temps numeric vector; Inverse temperature ladder, \eqn{t\in[0,1]}.
#' 
#' @param calibr.info list; Transformation parameters for
#' adjusting the pseudolikelihood function. See \code{\link[Bergm]{adjustPL}}. 
#'
#' @references
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
#' calibration.list <- list(Theta_MLE = info.adjustPL$Theta_MLE,
#'                          Theta_PL  = info.adjustPL$Theta_PL, 
#'                          W         = info.adjustPL$W, 
#'                          C         = info.adjustPL$C)
#'                          
#' # Specify location and shape of prior distribution: 
#' prior.mean  <- rep(0,2)
#' prior.sigma <- diag(5, 2)   
#'                                                 
#' # MCMC sampling and evidence estimation:
#' tempvec.powerp  <- seq(0, 1, length.out = 20)^5
#' 
#' pp.est.evidence <- evidence_PP(formula = flo.formula,
#'                                prior.mean = prior.mean,   
#'                                prior.sigma = prior.sigma,
#'                                nits = 10000,
#'                                burnin = 2000,
#'                                temps = tempvec.powerp,
#'                                calibr.info = calibration.list)
#'                                    
#' # MCMC diagnostics and posterior summaries:
#' bergm.output(pp.est.evidence)
#'   
#' # Log-marginal likelihood estimate:             
#' flo.model.logevidence <- pp.est.evidence$log.evidence
#'}
#'
#' @export
#'

evidence_PP <-function(formula,
                           prior.mean,   
                           prior.sigma,
                           nits,
                           burnin,
                           thin = 1,
                           tunePL = 2,
                           seed = 1,
                           temps = seq(0,1,length.out=50)^5,
                           calibr.info){
  
  #==========================
  # SUB-ROUTINES
  
  # 1. Obtain coefficient names for each model:
  mspecs <- function(x) {
      y <- ergm.getnetwork(x)
      model <- ergm_model(x, y)
      model$coef.names
  }
  
  # 2. Log pseudo-likelihood:
  expit <- function(x) exp(x) / (1 + exp(x)) 
  logPL.corr <- function(theta,
                         y,
                         X,
                         weights,
                         calibr.info){

    # Transform theta:
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL )
    xtheta <- c(X %*% theta_transf)
    #p <- expit(xtheta) # useless?
    log.like <- sum(dbinom(weights * y, weights, expit(xtheta), log = TRUE)) 
    
    return(log.like)
  }
  
  # 3. Log pseudo-posterior:
  logPP.corr <- function(theta,
                         y,
                         X,
                         weights,
                         prior.mean,   
                         prior.sigma,
                         temperature,
                         calibr.info){
    
    # Transform theta:
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL )
    
    # Calibrate only the pseudolikelihood:
    xtheta   <- c(X %*% theta_transf)
    p        <- expit(xtheta)
    log.like <- sum( dbinom(weights * y, weights, expit(xtheta), log = TRUE) )
    
    # Do not change the prior:
    log.prior<- dmvnorm(theta, mean = prior.mean, prior.sigma, log = TRUE) 
  
    return(temperature * log.like + log.prior)
    
  }
  
  score.temp.logpp <- function(theta,          
                               y,
                               X,
                               weights,
                               prior.mean,
                               prior.sigma,
                               calibr.info,
                               temperature){
    
    score.log.prior  <- -solve(prior.sigma, (theta - prior.mean))
    
    # Transform theta:
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL )
    
    p <- c( exp(as.matrix(X) %*% theta_transf ) / ( 1 + exp(as.matrix(X) %*% theta_transf ) ) )
    deriv.logpl <- c(t(X) %*% (weights*(y-p)))
    
    out <- temperature * (t(calibr.info$W) %*% deriv.logpl) + score.log.prior
    return(out)
    
  }
  
  #==========================
  # Hessian of corrected log pseudo-likelihood:
  Hessian.corrected.logpl <- function(theta,          
                                      y,
                                      X,
                                      weights,
                                      calibr.info){
    
    # Need to evaluate the Hessian of the psuedolikelihood at theta, not g(theta):
    p <- exp(as.matrix(X) %*% theta) / ( 1 + exp(as.matrix(X) %*% theta) )
    W <- Diagonal(x = as.vector( weights * p * (1 - p))) 
    
    Hessian <- -t(as.matrix(X)) %*% W %*% as.matrix(X)
    Hessian <- as.matrix(Hessian)
    
    return(t(calibr.info$W) %*% Hessian %*% calibr.info$W) 
  }
  
  #==========================
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
                        family  ="binomial"))
  
  # Model dimensions:
  dim <- length(mspecs(formula))
  
  # Tuning the variance-covariance matrix of the proposal distribution:
  tune.mat.PL <- diag(rep(tunePL, dim))
  B0 <- solve(prior.sigma)    
  
  H <- Hessian.corrected.logpl(calibr.info$Theta_PL,          
                               y = data.glm.initial[,"responses"],
                               X = data.glm.initial[,-c(1,2)],
                               weights = data.glm.initial[,"weights"],
                               calibr.info = calibr.info)
  
  PL.prop.sigma.mat <- tune.mat.PL %*% solve( B0 - H ) %*% tune.mat.PL 
  rownames(PL.prop.sigma.mat) <- colnames(PL.prop.sigma.mat) <-  rownames(Vcov.MPLE)
  
  #==========================
  # Output structures:
  numtemp <- length(temps)
  pplist <- list()
  acceptances <- rep(0,numtemp)
  cv.ell.D2 <- rep(0,numtemp)
  cv.vll.D2 <- rep(0,numtemp)
  htheta.D2 <- list()
  
  #==========================
  # Scale proposal variance within different temperatures:
  taup <- 1/sqrt( diag(PL.prop.sigma.mat) )
  tau0 <- solve(prior.sigma)
  
  a.tilde <- log(tau0[1,1] / taup) / log(temps[2]) 
  covar.temps <- list()
  
    for (i in numtemp:1) { 
      covar.temps[[i]]       <- PL.prop.sigma.mat
      # Take care of the variances:
      diag(covar.temps[[i]]) <- diag(PL.prop.sigma.mat) / (temps[i]^a.tilde[1])
    }
    
  # To avoid Inf variance at the zero temperature:
  covar.temps[[1]] <- covar.temps[[2]]
  
  #==========================
  # Define the progress bar for each temperature:
  #message("---MCMC start---")
  pb <- txtProgressBar(min = 0, max = numtemp, style = 3)  
  clock.start <- Sys.time() 
  
  #==========================
  # Run the sampler for each temperature:
  for(l in numtemp:1) { 
    if (l==1) {
      pplist[[l]]    <- rmvnorm(nits, mean = prior.mean, sigma = prior.sigma)
      acceptances[l] <- 1
    } else {
      if (l == numtemp) l0 <- calibr.info$Theta_PL else l0  <- c(tail(red)[6,])
      
      capture.output(red <- MCMCmetrop1R(logPP.corr,
                                         # MCMC parameters:
                                         mcmc       = nits, 
                                         burnin     = burnin,
                                         thin       = thin,
                                         theta.init = l0,
                                         seed       = seed,
                                         V          = covar.temps[[l]],
                                         # Model-related parameters:
                                         y           = mplesetup$response,
                                         X           = mplesetup$predictor,
                                         weights     = mplesetup$weights,
                                         calibr.info = calibr.info,
                                         prior.mean  = prior.mean,
                                         prior.sigma = prior.sigma,
                                         temperature = temps[l],
                                         logfun      = TRUE) )
      pplist[[l]] <- red
      acceptances[l] <- round(1 - rejectionRate(red)[1], 3)
    }
    setTxtProgressBar(pb, numtemp - l + 1) 
  }
  #==========================
  names(acceptances) <- paste("Temp", "_", 1:numtemp, " = ", round(temps, 5), sep = "")
  
  #==========================
  # Extract the mean and the 
  # variance of the log likelihoods:
  grad.logpp.list <- logPL.list <- list()
  
  #==========================
  # Define the progress bar for each temperature:
  #message("---CTI start---")
  pb <- txtProgressBar(min = 0, max = numtemp, style = 3)  
  
  for(j in 1:numtemp){
    
    # Need a special calculation at zero temperature:
    if (j == 1){
      grad.logpp.list[[j]]   <- t(apply(data.frame( pplist[[j]] ), 1,
                                        function(x) {
                                          -solve(prior.sigma, (x - prior.mean))
                                        } 
                                       )
                                  )
      
      logPL.list[[j]]    <- apply(pplist[[j]] , 1,
                                  function(x) {
                                    dmvnorm(x,
                                            prior.mean,
                                            prior.sigma,
                                            log = TRUE)
                                  } 
                                 )
    } else {
      grad.logpp.list[[j]] <- t(apply(data.frame( pplist[[j]] ), 1,
                                      function(x) {
                                        score.temp.logpp(theta       = x,          
                                                         y           = mplesetup$response,
                                                         X           = mplesetup$predictor,
                                                         weights     = mplesetup$weights,
                                                         prior.mean  = prior.mean,
                                                         prior.sigma = prior.sigma,
                                                         temperature = temps[j],
                                                         calibr.info = calibr.info)
                                                  }
                                      )
                                )
      logPL.list[[j]] <- apply(data.frame( pplist[[j]] ), 1,
                               function(x) {
                                 logPL.corr(theta      = x,          
                                            y          = mplesetup$response,
                                            X          = mplesetup$predictor,
                                            weights    = mplesetup$weights,
                                            calibr.info= calibr.info)
                                            }
                               ) 
    }
    
    ####################
    k <- dim * (dim + 3) / 2 
    l <- 2 * dim + 1         
    w.mat <- matrix(0, nrow = nits, ncol = k) 
    w.mat[, 1:dim] <- grad.logpp.list[[j]]         
    w.mat[, (dim + 1):(2 * dim)] <- as.matrix(pplist[[j]]) * grad.logpp.list[[j]] + rep(1, nits) 
    
    for (k1 in (1:(dim - 1)) ) {
      for (k2 in (k1 + 1):dim)  {
        w.mat[, l] <- as.matrix( pplist[[j]] )[, k1] * grad.logpp.list[[j]][, k2] + 
                     as.matrix( pplist[[j]] )[, k2] * grad.logpp.list[[j]][, k1]   
        l <- l + 1 
      }
    }

    phi.D2 <- c(-solve(var(w.mat)) %*% c(cov(w.mat, logPL.list[[j]])))
    htheta.D2[[j]] <- c(phi.D2 %*% t(w.mat))  
    
    cv.ell.D2[j] <- mean(logPL.list[[j]] + htheta.D2[[j]]) 
    cv.vll.D2[j] <- var(logPL.list[[j]] + htheta.D2[[j]]) 
    
    setTxtProgressBar(pb, j) 
    
  }
  clock.end<- Sys.time()
  runtime  <- difftime(clock.end, clock.start) 
  
  #==========================
  # Estimate the log marginal likelihood
  # using the standard and the modified
  # trapezoidal schemes:
  ppml <- function(cv.ell.D2, 
                   cv.vll.D2,
                   tempvec){
    
    N <- length(cv.ell.D2)
    cv.res.D2.mts <- 0

    for(i in 1:(N-1)){
      wts <- tempvec[i+1]-tempvec[i]

      # CTI-Modified trapezoidal rule, with quadratic polynomials:
      cv.res.D2.mts <- cv.res.D2.mts + wts*( (cv.ell.D2[i+1] + cv.ell.D2[i])/2.0 )-( (wts^2)/12 )*( cv.vll.D2[i+1] - cv.vll.D2[i] )
    }
    return( cv.res.D2.mts )
  }
  
  pp.estimates <- ppml(cv.ell.D2=cv.ell.D2, cv.vll.D2=cv.vll.D2, tempvec  =temps)
  log.evidence <- log(calibr.info$C) + pp.estimates
  
  #==========================
  ess <- round(effectiveSize(pplist[[numtemp]]),0)
  names(ess) <- mspecs(formula)
  
  #==========================
  out <- list(Theta = pplist[[numtemp]], # mcmc object; MCMC draws at temp_j=1 
              log.evidence = log.evidence, # CTI-Modified trapezoidal rule (quadratic polynomial)
              Time = runtime,
              formula = formula,
              AR_temp = acceptances,
              AR = acceptances[numtemp],
              ESS = ess,                 # Effective sample size for temperature t_j=1
              dim = dim,
              specs = mspecs(formula))
  
  class(out) <- "evidence_PP"
  return(out)
}






