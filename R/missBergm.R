#' Parameter estimation for Bayesian ERGMs under missing data
#'
#' Function to fit Bayesian exponential random graphs models under missing data
#' using the approximate exchange algorithm.
#'
#' @param formula formula; an \code{\link[ergm]{ergm}} formula object,
#' of the form  <network> ~ <model terms>
#' where <network> is a \code{\link[network]{network}} object
#' and <model terms> are \code{ergm-terms}.
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
#' @param prior.sigma square matrix; variance/covariance matrix for the multivariate Normal prior.
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
#' @param seed count;
#' random number seed for the Bergm estimation.
#'
#' @param startVals numeric matrix;
#' Starting values for the parameter estimation. startVals requires a matrix with parameters by number of chains. If nchains == NULL, nchains is equal to 2 * the number of parameters in the model.
#'
#' @param nImp count;
#' number of imputed networks to be returned. If null, no imputed network will be returned.
#'
#' @param missingUpdate count;
#' number of tie updates in each imputation step. By default equal to number of missing ties. Smaller numbers increase speed. Larger numbers lead to better sampling.
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
#' Koskinen, J.H., Robins, G.L., Pattison, P.E. (2010), "Analysing exponential
#' random graph (p-star) models with missing data using bayesian data augmentation,"
#' Statistical Methodology 7(3), 366-384
#'
#' Krause, R.W., Huisman, M., Steglich, C., Snijders, T.A. (2018), "Missing network
#' data a comparison of different imputation methods," Proceedings of the 2018
#' IEEE/ACM International Conference on Advances in Social Networks Analysis and
#' Mining 2018
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data(florentine)
#'
#'
#' # Create missing data
#' set.seed(22101992)
#'
#' missNode <- sample(1:16,1)
#' flomarriage[missNode,] <- NA
#' flomarriage[,missNode] <- NA
#'
#'
#' # Posterior parameter estimation:
#' m.flo <- missBergm(flomarriage ~ edges + kstar(2),
#'                burn.in = 50,
#'                aux.iters = 500,
#'                main.iters = 500,
#'                gamma = 1,
#'                nImp = 5)
#'
#' # Posterior summaries:
#'
#' bergm.output(m.flo)
#'
#' # Bayesian goodness-of-fit test:
#'
#' bgof(m.flo,
#'      aux.iters = 500,
#'      sample.size = 50,
#'      n.deg = 10,
#'      n.dist = 9,
#'      n.esp = 6)
#'}
#' @export
#'

missBergm <- function (formula,
                   burn.in = 100,
                   main.iters = 1000,
                   aux.iters = 1000,
                   prior.mean = NULL,
                   prior.sigma = NULL,
                   nchains = NULL,
                   gamma = 0.5,
                   sigma.epsilon = NULL,
                   seed = NULL,
                   startVals = NULL,
                   nImp = NULL,
                   missingUpdate = NULL,
                   ...) {

    if (is.null(seed)) {
      set.seed(sample(1:999999,1))
    } else {
      set.seed(seed)
    }

    y <- ergm.getnetwork(formula)
    model <- ergm_model(formula, y)
    sy <- summary(formula)
    dim <- length(sy)

    # ---
    impNets <-  NULL

    if (!any(is.na(as.matrix.network(y)))) {
      stop("Network has no missing data. Use bergm() instead.")
    }

    if (!is.null(nImp)) {
      nImp <- max(0,min(nImp, main.iters))
      thinImp <- as.integer((main.iters)/nImp)
      impIter <- 1
      impNets <- vector("list",nImp)
    }

    missingTies <- matrix(0,y$gal$n,y$gal$n)
    missingTies[is.na(as.matrix.network(y))] <- 1
    missingTies <- as.edgelist(as.network(missingTies), n = y$gal$n)

    if (is.null(missingUpdate)) {missingUpdate <- sum(is.na(as.matrix.network(y)))}


    # ---
    Clist <- ergm.Cprepare(y, model)

    control <- control.ergm(MCMC.burnin = aux.iters,
                            MCMC.interval = 1,
                            MCMC.samplesize = 1)

    proposal <- ergm_proposal(object = ~.,
                              constraints = ~.,
                              arguments = control$MCMC.prop.args,
                              nw = y)
    # ---

    snooker <- 0
    if (is.null(prior.mean)) prior.mean <- rep(0, dim)
    if (is.null(prior.sigma)) prior.sigma <- diag(100, dim)
    if (is.null(nchains)) nchains <- 2 * dim
    if (is.null(sigma.epsilon)) sigma.epsilon <- diag(0.0025, dim)
    if (dim == 1) {
    	 nchains <- 1
       sigma.epsilon <- diag(gamma, dim)
    }
    Theta <- array(NA, c(main.iters, dim, nchains))


    if (is.null(startVals)) {
      theta <- matrix(runif(dim * nchains, min = -0.1, max = 0.1), dim, nchains)
    } else if (nrow(startVals) != Clist$nstats || ncol(startVals) != nchains) {
      stop("StartVals has wrong dimensions. Startvals requires a matrix with nrow = nmber of parameters and ncol = number of chains. If nchains == NULL, nchains is 2 * number of parameters.")
    } else {
      theta <- startVals
    }

    acc.counts <- rep(0L, nchains)
    theta1 <- rep(NA, dim)
    tot.iters <- burn.in + main.iters


    impNet <- y

    f <- as.character(formula)
    currentFormula <- formula(paste("impNet",
                                    f[3:length(f)], sep = " ~ "))

    clock.start <- Sys.time()

    for (k in 1:tot.iters) {
        for (h in 1:nchains) {

            # --- Step 1 Parameter proposal
            if (dim > 1 && nchains > 1) {
                snooker <- gamma * apply(theta[, sample(seq(1, nchains)[-h], 2)], 1, diff)
            }
            theta1 <- theta[, h] + snooker + rmvnorm(1, sigma = sigma.epsilon)[1, ]

            # --- Step 2 Network simulation
            delta <- ergm_MCMC_slave(Clist = Clist,
                                     proposal = proposal,
                                     eta = theta1,
                                     control = control,
                                     verbose = FALSE)$s
            # --- Step 3 Proposal acceptence/rejection

            pr <- dmvnorm(rbind(theta1, theta[, h]),
                          mean = prior.mean,
                          sigma = prior.sigma,
                          log = TRUE)

            beta <- (theta[, h] - theta1) %*% t(delta) + pr[1] - pr[2]

            if (beta >= log(runif(1))) {
              theta[, h] <- theta1

              if (k > burn.in) { acc.counts[h] <- acc.counts[h] + 1 }

              # --- Step 4 Network Imputation
              impNet <- simulate(currentFormula, coef = theta1,
                                           statsonly = FALSE,
                                           basis = y, constraints =
                                             ~fixallbut(missingTies),
                                           nsim = 1,
                                           control = control.simulate(
                                             MCMC.burnin = missingUpdate))

              y2 <- ergm.getnetwork(currentFormula)
              model2 <- ergm_model(currentFormula, y2)
              Clist <- ergm.Cprepare(y2, model2)
            }


          }

        if (k > burn.in) Theta[k - burn.in, , ] <- theta

         if (!is.null(nImp)) {
          if ((k - burn.in) == impIter * thinImp) {
            impNets[[impIter]] <- impNet
            impIter <- impIter + 1
          }
        }
    }



    if (nchains == 1) Theta <- as.matrix(Theta[, , 1])
    if (nchains > 1) Theta <- apply(Theta, 2, cbind)

    clock.end <- Sys.time()
    runtime <- difftime(clock.end, clock.start)

    out = list(Time = runtime,
               formula = formula,
               model = model,
               specs = model$coef.names,
               dim = dim,
               nchains = nchains,
               stats = sy,
               Theta = Theta,
               AR = acc.counts / main.iters,
               impNets = impNets)
    class(out) <- "bergm"
    return(out)
  }