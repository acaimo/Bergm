#' Adjustment of ERGM pseudolikelihood
#'
#' Function to estimate the transformation parameters for
#' adjusting the pseudolikelihood function.
#' 
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}.
#' 
#' @param aux.iters count; number of auxiliary iterations used for drawing the first network from the ERGM likelihood. See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param n.aux.draws count; Number of auxiliary networks drawn from the ERGM likelihood. See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param aux.thin count; Number of auxiliary iterations between network draws after the first network is drawn. See \code{\link[ergm]{control.simulate.formula}}.
#' 
#' @param ladder count; Length of temperature ladder (>=3).
#' 
#' @param estimate If "MLE" (the default), then an approximate maximum likelihood estimator is returned. If "CD" , the Monte-Carlo contrastive divergence estimate is returned. See \code{\link[ergm]{ergm}}.
#' 
#' @param seed integer; seed for the random number generator. See \code{set.seed}.
#' 
#' @param ... Additional arguments, to be passed to the ergm function. See \code{\link[ergm]{ergm}}.
#'
#' @references
#' Bouranis, L., Friel, N., & Maire, F. (2018). Bayesian model selection for exponential 
#' random graph models via adjusted pseudolikelihoods. 
#' Journal of Computational and Graphical Statistics, 27(3), 516-528. \url{https://arxiv.org/abs/1706.06344}
#'
#' @export
#'
ergmAPL <- function(formula, 
                    aux.iters   = NULL, 
                    n.aux.draws = NULL, 
                    aux.thin    = NULL, 
                    ladder      = NULL,
                    estimate    = c("MLE","CD"),
                    seed        = 1,
                    ...) 
{
  y     <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  n     <- dim(y[,])[1]
  sy    <- summary(formula)
  dim   <- length(sy)
  
  if (dim == 1) stop("Model dimension must be greater than 1")
  
  if (is.null(aux.iters)) aux.iters <- 3000
  if (is.null(n.aux.draws)) n.aux.draws <- 50
  if (is.null(aux.thin)) aux.thin <- 50
  if (is.null(ladder)) ladder <- 50
  
  estimate <- match.arg(estimate)
  
  control <- control.ergm(MCMC.burnin        = aux.iters, 
                          MCMC.interval      = aux.thin, 
                          MCMC.samplesize    = n.aux.draws, 
                          MCMC.maxedges      = Inf,
                          seed               = seed)
  
  adjusted_logPL <- function(theta, Y, X, weights, calibr.info) {
    theta_transf   <- c(calibr.info$W %*% (theta - calibr.info$Theta_MLE) + calibr.info$Theta_PL)
    xtheta         <- c(X %*% theta_transf)
    log.like       <- sum(dbinom(weights * Y, weights, exp(xtheta)/(1 + exp(xtheta)), log = TRUE))
    return(log.like)
  }
  
  Hessian_logPL <- function(theta, X, weights) {
    p <- exp(as.matrix(X) %*% theta)/(1 + exp(as.matrix(X) %*% theta))
    W <- Diagonal(x = as.vector(weights * p * (1 - p)))
    Hessian <- -t(as.matrix(X)) %*% W %*% as.matrix(X)
    Hessian <- as.matrix(Hessian)
    return(Hessian)
  }
  
  mplesetup <- ergmMPLE(formula)
  data.glm.initial <- cbind(mplesetup$response, 
                            mplesetup$weights, 
                            mplesetup$predictor)
  colnames(data.glm.initial) <- c("responses", "weights", 
                                  colnames(mplesetup$predictor))
  path <- seq(0, 1, length.out = ladder)
  
  suppressMessages(mle  <- ergm(formula, estimate = estimate, verbose = FALSE, control = control.ergm(seed = seed), ...))
  suppressMessages(mple <- ergm(formula, estimate = "MPLE", verbose = FALSE))
  
  if( any( c(-Inf, Inf) %in% mle$coefficients) | 
      any( c(-Inf, Inf) %in% mple$coefficients) ) {
    
    if( any( c(-Inf, Inf) %in% mle$coefficients) ){
      
      inf_term <- names(mle$coefficients[mle$coefficients %in% c(-Inf, Inf)])
      inf_term <- paste(inf_term, collapse = " & ")
      
    } else if( any( c(-Inf, Inf) %in% mple$coefficients) ){
    
      inf_term <- names(mple$coefficients[mple$coefficients %in% c(-Inf, Inf)])
      inf_term <- paste(inf_term, collapse = " & ")
    }
    
    stop( paste0("The MLE/MPLE contains an infinite value for the following model terms: ", 
                 inf_term, 
                 ". Consider changing these model terms.") )
  }

  y0 <- simulate(formula, 
                 coef        = mle$coefficients, 
                 nsim        = 1, 
                 control     = control.simulate(MCMC.burnin = 1,
                                                MCMC.interval = 1),
                 return.args = "ergm_state")$object
  
  z <- as.matrix( ergm_MCMC_sample(y0,
                                   theta    = mle$coefficients,
                                   stats0   = sy,
                                   control  = control
                                   )$stats[[1]] )
  
  H   <- -cov(z)
  
  HPL <- Hessian_logPL(theta   = mple$coefficients, 
                       X       = mplesetup$predictor, 
                       weights = mplesetup$weights)

  HPL <- round(HPL,5)
  mat <- -H           
  mat <- round(mat,5)
  
  if( is.positive.definite(mat) == TRUE & is.positive.definite(HPL) == TRUE ) {
    chol.HPL<- chol(-HPL)
    chol.H  <- chol(mat)
    W       <- solve(chol.HPL) %*% chol.H
    
  } else if (is.positive.definite(mat) == FALSE & is.positive.definite(HPL) == TRUE){
    
    suppressWarnings( chol.H <- chol(mat, pivot = TRUE) )
    pivot    <- attr(chol.H, "pivot")
    oo       <- order(pivot)
    chol.HPL <- chol(-HPL)
    W        <- solve(chol.HPL) %*% chol.H[, oo]
    
  } else if (is.positive.definite(mat) == TRUE & is.positive.definite(HPL) == FALSE){
    chol.H <- chol(mat)
    
    suppressWarnings( chol.HPL <- chol(-HPL, pivot = TRUE) )
    pivot <- attr(chol.HPL, "pivot")
    oo    <- order(pivot)
    W     <- solve(chol.HPL[, oo]) %*% chol.H
    
  } else {
    suppressWarnings( chol.H <- chol(mat, pivot = TRUE) )
    pivot.H <- attr(chol.H, "pivot")
    oo.H    <- order(pivot.H)
    
    suppressWarnings( chol.HPL <- chol(-HPL, pivot = TRUE) )
    pivot.HPL <- attr(chol.HPL, "pivot")
    oo.HPL    <- order(pivot.HPL)
    W         <- solve(chol.HPL[, oo.HPL]) %*% chol.H[, oo.H]
  }
  
  adjust.info <- list(Theta_MLE = mle$coefficients, 
                      Theta_PL  = mple$coefficients, 
                      W = W)
  
  E <- lapply(seq_along(path)[-length(path)], function(i) { 
    
    y0E <- simulate(formula, 
                    coef        = path[i] * adjust.info$Theta_MLE, 
                    nsim        = 1, 
                    control     = control.simulate(MCMC.burnin   = 1,
                                                   MCMC.interval = 1),
                   return.args = "ergm_state")$object
    
    Monte.Carlo.samples <- as.matrix( ergm_MCMC_sample(y0E,
                                                       theta    = path[i] * adjust.info$Theta_MLE,
                                                       stats0   = sy,
                                                       control  = control)$stats[[1]])
    
    log( mean( exp( (path[i + 1] - path[i]) * adjust.info$Theta_MLE %*% t(Monte.Carlo.samples) ) ) )
  })
  
  E <- sum(unlist(E))
  
  if (y$gal$directed == FALSE) logz0A <- choose(n, 2) * log(2) else logz0A <- (n * (n - 1)) * log(2)
  
  logztheta <- E + logz0A
  
  ll.true <- c(matrix(adjust.info$Theta_MLE, nrow = 1) %*% sy) - logztheta
  
  ll.adjpseudo <- adjusted_logPL(theta       = adjust.info$Theta_MLE, 
                                 Y           = mplesetup$response, 
                                 X           = mplesetup$predictor, 
                                 weights     = mplesetup$weights, 
                                 calibr.info = adjust.info)
  
  out <- list(formula      = formula, 
              Theta_MLE    = mle$coefficients, 
              Theta_PL     = mple$coefficients, 
              W            = W, 
              logC         = ll.true - ll.adjpseudo, 
              ll_true      = ll.true, 
              logztheta    = logztheta, 
              ll_adjpseudo = ll.adjpseudo)
  return(out)
}