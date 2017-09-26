#' Bayesian parameter inference for ERGMs
#'
#' Function to fit Bayesian exponential random graphs models
#' using the approximate exchange algorithm.
#'
#' @param formula formula; an \code{R} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link{network}} object
#' and <model terms> are \code{\link{ergm-terms}}.
#'
#' @param burn.in count; number of burn-in iterations at the beginning of an MCMC run.
#' If \code{nchains} > 2,
#' it refers to the number of burn-in iterations for every chain of the population.
#'
#' @param main.iters count; number of iterations for the MCMC chain(s) excluding burn-in.
#' If \code{nchains} > 2,
#' it refers to the number of iterations for every chain of the population.
#'
#' @param aux.iters count; number of auxiliary iterations used for network simulation.
#'
#' @param prior.mean vector; mean vector of the multivariate Normal prior.
#' By default set to a vector of 0's.
#'
#' @param sigma.mean square matrix; variance/covariance matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100.
#'
#' @param nchains count; number of chains of the population MCMC.
#' By default set to twice the model dimension (number of model terms).
#' If the model is one-dimensional, \code{nchains} = 1.
#'
#' @param gamma scalar; ``parallel ADS move factor''.
#' If the model is one-dimensional, \code{nchains} = 1 and
#' \code{gamma} = \code{sigma.espilon}
#' and is used as the variance of the Normal proposal distribution.
#'
#' @param sigma.epsilon square matrix;
#' variance/covariance matrix for the multivariate Normal proposal when \code{nchains} > 2.
#' By default set to a diagonal matrix with every diagonal entry equal to 0.0025.
#' If the model is one-dimensional, \code{sigma.espilon} = \code{gamma}
#' and is used as the variance of the Normal proposal distribution.
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
#'
#' p.flo <- bergm(flomarriage ~ edges + kstar(2),
#'                burn.in = 50,
#'                aux.iters = 500,
#'                main.iters = 500,
#'                gamma = 1)
#'
#' # MCMC diagnostics and posterior summaries:
#'
#' bergm.output(p.flo)
#'
#' # Bayesian goodness-of-fit test:
#'
#' bgof(p.flo,
#'      aux.iters = 500,
#'      sample.size = 50,
#'      n.deg = 10,
#'      n.dist = 9,
#'      n.esp = 6)
#'
#' @import network
#' @import ergm
#' @import mvtnorm
#'
#' @export
#'

bergm <- function (formula,
                   burn.in = 100,
                   main.iters = 1000,
                   aux.iters = 1000,
                   prior.mean = NULL,
                   sigma.mean = NULL,
                   nchains = NULL,
                   gamma = 0.5,
                   sigma.epsilon = NULL,
                   
                   ...){

    y <- ergm.getnetwork(formula)
    model <- ergm.getmodel(formula, y)
    Clist <- ergm.Cprepare(y, model)

    stats0 <- summary(formula)
    control <- control.simulate.formula(MCMC.burnin = aux.iters,
                                        MCMC.interval = 0)
    control$MCMC.samplesize <- 1

    MHproposal <- MHproposal.ergm(object = model,
                                  constraints = ~., arguments = control$MCMC.prop.args,
                                  nw = y, weights = control$MCMC.prop.weights, class = "c",
                                  reference = ~Bernoulli, response = NULL)

    snooker <- 0
    if (is.null(prior.mean)) prior.mean <- rep(0, Clist$nstats)
    if (is.null(sigma.mean)) sigma.mean <- diag(100, Clist$nstats)
    if (is.null(nchains)) nchains <- 2 * Clist$nstats
    if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, Clist$nstats)
    if (Clist$nstats == 1) {
    	 nchains <- 1
        sigma.epsilon <- diag(gamma, Clist$nstats)
    }
    Theta <- array(NA, c(main.iters, Clist$nstats, nchains))
    theta <- matrix(runif(Clist$nstats * nchains, min = -0.1, max = 0.1), Clist$nstats, nchains)
    acc.counts <- rep(0L, nchains)
    theta1 <- rep(NA, Clist$nstats)
    tot.iters <- burn.in + main.iters

    for (k in 1L:tot.iters) {
        for (h in 1L:nchains) {
            if (Clist$nstats > 1 && nchains > 1) {
                snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
            }

            theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1, ]

            pr <- dmvnorm(rbind(theta1, theta[, h]),
                          mean = prior.mean,
                          sigma = sigma.mean,
                          log=TRUE)

            delta <- ergm.mcmcslave(Clist,
                                    MHproposal,
                                    eta0 = theta1,
                                    control,
                                    verbose = FALSE)$s

            beta <- (theta[, h] - theta1) %*% t(delta) + pr[1] - pr[2]

            if (beta >= log(runif(1))) {

                theta[, h] <- theta1
                if (k > burn.in) acc.counts[h] <- acc.counts[h] + 1

            }

        }
        if (k > burn.in) Theta[k - burn.in, , ] <- theta
    }
    if (nchains == 1) Theta <- as.matrix(Theta[, , 1])
    if (nchains > 1) Theta <- apply(Theta, 2, cbind)
    
    out = list(Clist = Clist, MHproposal = MHproposal, control = control,
        formula = formula, model = model, nnodes = Clist$n, specs = model$coef.names,
        dim = Clist$nstats, nchains = nchains, stats = stats0,
        Theta = Theta, nchains = nchains, AR = acc.counts / main.iters,
        prior.mean = prior.mean, sigma.mean = sigma.mean, aux.iters = aux.iters)
    class(out) <- "bergm"
    return(out)

}
