#' Adjustment of pseudolikelihood function
#'
#' Function to estimate the transformation parameters for
#' adjusting the pseudolikelihood function.
#' 
#' @param ergm.formula formula; an \code{R} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link{network}} object
#' and <model terms> are \code{\link{ergm-terms}}.
#' 
#' @param aux.iters count; Number of proposals before any MCMC sampling is done. See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param noisy.nsim count; Number of TNT draws. See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param noisy.thin count; Number of proposals between sampled statistics. See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param ladder count; Length of temperature ladder (>=3).
#' 
#' @param ... Additional arguments, to be passed to the ergm function. See \code{\link[ergm]{ergm}}.
#'
#' @references
#' Bouranis, L., Friel, N., and Maire, F. (2017). Bayesian model selection for exponential random graph models via
#' adjusted pseudolikelihoods. \url{https://arxiv.org/abs/1706.06344}
#'
#' @export
#'

adjustPL <- function(ergm.formula,
                     aux.iters    = 3000,
                     noisy.nsim   = 50,  
                     noisy.thin   = 50,  
                     ladder       = 50,  
                     ...                  
                     ){

  #==========================
  y <- ergm.getnetwork(ergm.formula) 

  #==========================
  model.sy <- summary(ergm.formula)

  #==========================
  expit <- function(x) exp(x) / (1 + exp(x)) 
  
  #==========================
  logPL.corr <- function(theta,
                         y,
                         X,
                         weights,
                         calibr.info){
    
    theta_transf <- c( calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL ) 
    
    xtheta   <- c(X%*%theta_transf)
    p        <- expit(xtheta)
    log.like <- sum( dbinom(weights*y, weights, expit(xtheta), log = TRUE) )
    
    return(log.like)
    
  }# End function
  
  #==========================
  # Hessian of log pseudolikelihood:
  Hessian.logPL <-function(theta,          
                           y,
                           X,
                           weights)
  {
    
    p <- exp(as.matrix(X)%*%theta)/ ( 1+exp(as.matrix(X)%*%theta) )
    W <- Diagonal( x=as.vector( weights*p*(1-p) )) 
    
    Hessian <- - t( as.matrix(X) ) %*% W %*% as.matrix(X)
    Hessian <- as.matrix(Hessian)
    
    return( Hessian ) 
  }# End function
  
  #==========================
  # Get data in aggregated format:
  mplesetup                 <- ergmMPLE(ergm.formula)
  data.glm.initial          <- cbind(mplesetup$response, 
                                     mplesetup$weights, 
                                     mplesetup$predictor)
  colnames(data.glm.initial) <- c("responses", 
                                  "weights", 
                                  colnames(mplesetup$predictor))
  
  # Variance-covariance matrix at MPLE:
  Vcov.MPLE <- vcov(glm(mplesetup$response ~ . - 1, 
                        data    = data.frame(mplesetup$predictor), 
                        weights = mplesetup$weights, 
                        family  = "binomial"))

  #==========================
  path <- seq(0, 1, length.out = ladder)
  
  #==========================
  post   <- lapply(list(ergm.formula), 
                   FUN=function(x) {
                     y       <- ergm.getnetwork(x)
                     model   <- ergm.getmodel(x, y)
                     Clist   <- ergm.Cprepare(y, model)
                     control <- control.simulate.ergm(MCMC.burnin       = aux.iters, 
                                                      MCMC.interval     = noisy.thin,
                                                      MCMC.init.maxedges= 20000,
                                                      MCMC.max.maxedges = Inf
                     )
                     
                     control$MCMC.samplesize <- noisy.nsim 
                     
                     MHproposal <- MHproposal.ergm(object      = model, 
                                                   constraints = ~., 
                                                   arguments   = control$MCMC.prop.args, 
                                                   nw          = y, 
                                                   weights     = control$MCMC.prop.weights, 
                                                   class       = "c",
                                                   reference   = ~Bernoulli, 
                                                   response    = NULL)
                     
                     list(Clist      = Clist, 
                          MHproposal = MHproposal, 
                          control    = control
                     )
                   }# End function
  )  
  
  #==========================
  # Model dimensions:
  mdims <- post[[1]]$Clist$nstats
  
  #==========================
  n <- post[[1]]$Clist$n
  
  #==========================
  message("---Mode estimation---")
  
  #==========================
  mle <- ergm(ergm.formula,...)
  
  #==========================
  # Obtain the MPLE for the model of interest:
  capture.output(mple <- ergm(ergm.formula, estimate= "MPLE") )
  
  message("---Curvature Adjustment---")
  
  #==========================
  # Simulate graphs at the MLE:
  sim.samples<- simulate(ergm.formula, 
                         coef     = mle$coef,
                         statsonly= T,
                         nsim     = noisy.nsim,
                         control  = control.simulate(MCMC.burnin   = aux.iters, 
                                                     MCMC.interval = noisy.thin) )
  
  #==========================
  # Hessian of true log-posterior: 
  Hessian.true.logLL <- -cov(sim.samples)
  
  HPL <- Hessian.logPL(theta   = mple$coef,          
                       y       = mplesetup$response,
                       X       = mplesetup$predictor,
                       weights = mplesetup$weights)
  
  #==========================
  chol.true.Hessian<- chol(- Hessian.true.logLL)
  chol.PL.Hessian  <- chol(- HPL) 
  
  #==========================
  # Calculate transformation matrix W:
  W <- solve(chol.PL.Hessian)%*%chol.true.Hessian

  adjust.info <- list(Theta_MLE= mle$coef,
                      Theta_PL = mple$coef,
                      W        = W)
  
  #==========================
  message("---Estimating log normalising constant at the MLE---")
  
  E<- lapply(seq_along(path)[-length(path)], 
               function(i){
                 Monte.Carlo.samples   <- ergm.mcmcslave(Clist      = post[[1]]$Clist, 
                                                         MHproposal = post[[1]]$MHproposal, 
                                                         eta0       = path[i]*adjust.info$Theta_MLE, 
                                                         control    = post[[1]]$control, 
                                                         verbose    = FALSE)$s
                 Monte.Carlo.samples   <- sweep(Monte.Carlo.samples, 2, model.sy, `+`)
                 
                 log( mean( exp( (path[i+1]-path[i])*adjust.info$Theta_MLE %*%t(Monte.Carlo.samples) ) ) )
                 
               })
  
  E <- sum(unlist(E))
  
  #==========================
  # Estimate log z(0) depending on the type of network:
  if (y$gal$directed==FALSE) logz0 <- choose(n,2)*log(2) else logz0 <- (n*(n-1))*log(2)
  
  logztheta <- E + logz0 #E + choose(n,2)*log(2)
               
  #==========================
  ll.true <- c( matrix(adjust.info$Theta_MLE, nrow=1)%*%model.sy ) - logztheta
  
  ll.adjpseudo <-  logPL.corr(theta      = adjust.info$Theta_MLE,
                              y          = mplesetup$response,
                              X          = mplesetup$predictor,
                              weights    = mplesetup$weights,
                              calibr.info= adjust.info)  
  
  C <- exp(ll.true - ll.adjpseudo )
  
  #==========================
  out<- list(Theta_MLE    = mle$coef,
             Theta_PL     = mple$coef,
             W            = W,
             C            = C,
             ll_true      = ll.true,
             logztheta    = logztheta,
             ll_adjpseudo = ll.adjpseudo)
  
  class(out) <-"adjustPL"
  
  return(out)

}# End function