#' Parameter estimation for Bayesian ERGMs
#'
#' Function to fit Bayesian exponential random graphs models
#' using the approximate exchange algorithm.
#'
#' @param formula formula; 
#' an \code{\link[ergm]{ergm}} formula object,
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
#' @param burn.in count; 
#' number of burn-in iterations for every chain of the population.
#'
#' @param main.iters count; 
#' number of iterations for every chain of the population.
#'
#' @param aux.iters count; 
#' number of auxiliary iterations used for network simulation.
#'
#' @param nchains count; 
#' number of chains of the population MCMC.
#' By default set to twice the model dimension (number of model terms).
#'
#' @param gamma scalar; 
#' parallel adaptive direction sampling move factor.
#' 
#' @param V.proposal count; 
#' diagonal entry for the multivariate Normal proposal.
#' By default set to 0.0025.
#'
#' @param startVals vector;
#' optional starting values for the parameter estimation. 
#'
#' @param ... additional arguments, to be passed to lower-level functions.
#'
#' @references
#' Caimo, A. and Friel, N. (2011), "Bayesian Inference for Exponential Random Graph Models,"
#' Social Networks, 33(1), 41-55. \url{http://arxiv.org/abs/1007.5192}
#'
#' Caimo, A. and Friel, N. (2014), "Bergm: Bayesian Exponential Random Graphs in R,"
#' Journal of Statistical Software, 61(2), 1-25. \url{jstatsoft.org/v61/i02}
#'
#' @examples
#' # Load the florentine marriage network
#' data(florentine)
#'
#' # Posterior parameter estimation:
#' p.flo <- bergm(flomarriage ~ edges + kstar(2),
#'                burn.in    = 50,
#'                aux.iters  = 500,
#'                main.iters = 1000,
#'                gamma      = 1.2)
#'
#' # Posterior summaries:
#' summary(p.flo)
#'
#' @export

bergm <- function(formula,
                  prior.mean  = NULL,
                  prior.sigma = NULL,
                  burn.in     = 100,
                  main.iters  = 1000,
                  aux.iters   = 1000,
                  nchains     = NULL,
                  gamma       = 0.5,
                  V.proposal  = 0.0025,
                  startVals   = NULL,
                  ...){

  y     <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  sy    <- summary(formula)
  dim   <- length(sy)
  
  if (dim == 1) stop("Model dimension must be greater than 1")
  if (any(is.na(as.matrix.network(y)))) print("Network has missing data. Use bermM() instead.")
  
  # For network simulation
  Clist   <- ergm.Cprepare(y, model)
  
  control <- control.ergm(MCMC.burnin = aux.iters,
                          MCMC.interval = 1,
                          MCMC.samplesize = 1)
  
  proposal <- ergm_proposal(object = ~.,
                            constraints = ~.,
                            arguments = control$MCMC.prop.args,
                            nw = y)

  if (is.null(prior.mean))  prior.mean  <- rep(0, dim)
  if (is.null(prior.sigma)) prior.sigma <- diag(100, dim, dim)
  if (is.null(nchains)) nchains <- 2 * dim
  
  S.prop <- diag(V.proposal, dim, dim)
  
  Theta <- array(NA, c(main.iters, dim, nchains))

  # starting values
  if (is.null(startVals)) {
    suppressMessages( mple <- ergm(formula, estimate = "MPLE", verbose = FALSE)$coef )
    theta <- matrix( mple + runif(dim * nchains, min = -0.1, max = 0.1), dim, nchains)
  } else {
    theta <- matrix( startVals + runif(dim * nchains, min = -0.1, max = 0.1), dim, nchains)
  }

  theta1     <- rep(NA, dim)
  tot.iters  <- burn.in + main.iters
  
  clock.start <- Sys.time()

  message(" > MCMC start")
  
  for (k in 1:tot.iters) {
      for (h in 1:nchains) {

          theta1 <- theta[, h] + 
                    gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff) + 
                    rmvnorm(1, sigma = S.prop)[1, ]
          
          delta <- ergm_MCMC_slave(Clist = Clist,
                                   proposal = proposal,
                                   eta = theta1,
                                   control = control,
                                   verbose = FALSE)$s
          
          pr <- dmvnorm(rbind(theta1, theta[, h]),
                        mean = prior.mean,
                        sigma = prior.sigma,
                        log = TRUE)
            
          beta <- (theta[, h] - theta1) %*% t(delta) + pr[1] - pr[2]

          if (beta >= log(runif(1))) theta[, h] <- theta1
      }
      if (k > burn.in) Theta[k - burn.in, , ] <- theta
  }

  clock.end <- Sys.time()
  runtime   <- difftime(clock.end, clock.start)
  
  FF <- mcmc(apply(Theta, 2, cbind))
  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL
  ess <- round(effectiveSize(FF), 0)
  names(ess) <- model$coef.names
  
  out = list(Time    = runtime,
             formula = formula,
             specs   = model$coef.names,
             dim     = dim,
             Theta   = FF,
             AR      = AR,
             ess     = ess)

  class(out) <- "bergm"
  return(out)
}
