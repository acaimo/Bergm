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
#' @param burn.in count; number of burn-in iterations for every chain of the population.
#'
#' @param main.iters count; number of iterations for every chain of the population.
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
#'
#' @param gamma scalar; parallel adaptive direction sampling move factor.
#'
#' @param V.proposal count; diagonal entry for the multivariate Normal proposal.
#' By default set to 0.0025.
#'
#' @param seed count;
#' random number seed for the Bergm estimation.
#'
#' @param startVals vector;
#' optional starting values for the parameter estimation. 
#'
#' @param offset.coef vector;
#' A vector of coefficients for the offset terms.
#' 
#' @param nImp count;
#' number of imputed networks to be returned. If null, no imputed network will be returned.
#'
#' @param missingUpdate count;
#' number of tie updates in each imputation step. 
#' By default equal to number of missing ties. 
#' Smaller numbers increase speed. Larger numbers lead to better sampling.
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
#' Statistical Methodology 7(3), 366-384.
#'
#' Krause, R.W., Huisman, M., Steglich, C., Snijders, T.A. (2018), "Missing network
#' data a comparison of different imputation methods," Proceedings of the 2018
#' IEEE/ACM International Conference on Advances in Social Networks Analysis and
#' Mining 2018.
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data(florentine)
#' 
#' # Create missing data
#' set.seed(22101992)
#' n <- dim(flomarriage[, ])[1]
#' missNode <- sample(1:n, 1)
#' flomarriage[missNode, ] <- NA
#' flomarriage[, missNode] <- NA
#'
#' # Posterior parameter estimation:
#' m.flo <- bergmM(flomarriage ~ edges + kstar(2),
#'                 burn.in    = 50,
#'                 aux.iters  = 500,
#'                 main.iters = 1000,
#'                 gamma      = 1.2,
#'                 nImp       = 5)
#'
#' # Posterior summaries:
#' summary(m.flo)
#'}
#' @export

bergmM <- function(formula, 
                   burn.in = 100, 
                   main.iters = 1000, 
                   aux.iters = 1000,
                   prior.mean = NULL, 
                   prior.sigma = NULL, 
                   nchains = NULL, 
                   gamma = 0.5,
                   V.proposal = 0.0025, 
                   seed = NULL, 
                   startVals = NULL,  
                   offset.coef = NULL,
                   nImp = NULL,
                   missingUpdate = NULL,
                   ...)
{
  if (is.null(seed))
    set.seed(sample(1:999, 1))
  else set.seed(seed)
  y <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  sy <- summary(formula)
  dim <- length(sy)
  if (dim == 1) {stop("Model dimension must be greater than 1")}
  if (!is.null(offset.coef)) {
    if (any(offset.coef %in% c(Inf,-Inf,NaN,NA))) {
      stop("Inf,-Inf,NaN,NA are not allowed for offset.coef. \n If Inf or -Inf are required use large values instead (e.g., 1000 or -1000).")
    }
  }
  if (!all(is.na(as.matrix.network(y)))) {
        print("Network has no missing data. Use bergm() for faster estimation instead.")}
  impNets <- NULL
  if (!is.null(nImp)) {
    nImp <- max(0, min(nImp, main.iters))
    thinImp <- as.integer(main.iters/nImp)
    impIter <- 1
    impNets <- vector("list", nImp)
  }
  missingTies <- matrix(0, y$gal$n, y$gal$n)
  missingTies[is.na(as.matrix.network(y))] <- 1
  missingTies <- as.edgelist(as.network(missingTies), n = y$gal$n)
  if (is.null(missingUpdate)) {
    missingUpdate <- sum(is.na(as.matrix.network(y)))
  }
  Clist <- ergm.Cprepare(y, model)
  control <- control.ergm(MCMC.burnin = aux.iters, MCMC.interval = 1,
                          MCMC.samplesize = 1)

  if (!is.null(control$init)) {
    if (length(control$init) != length(model$etamap$offsettheta)) {
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters.", "If you are passing output from another ergm run as control$init,",
                 "in a model with curved terms, see help(enformulate.curved)."))
    }
  } else {control$init <- rep(NA, length(model$etamap$offsettheta))}
  if (!is.null(offset.coef)) {
    if (length(control$init[model$etamap$offsettheta]) !=
        length(offset.coef)) {
      stop("Invalid offset parameter vector offset.coef: ",
           "wrong number of parameters: expected ",
           length(control$init[model$etamap$offsettheta]),
           " got ", length(offset.coef), ".")}
    control$init[model$etamap$offsettheta] <- offset.coef
  }
  if (any(is.na(control$init) & model$etamap$offsettheta)) {
    stop("The model contains offset terms whose parameter values have not been specified:",
         paste.and(model$coef.names[is.na(control$init) |
                                model$offsettheta]), ".", sep = "") }


  proposal <- ergm_proposal(object = ~., constraints = ~.,
                            arguments = control$MCMC.prop.args, nw = y)
  if (is.null(prior.mean))
    prior.mean <- rep(0, dim)
  if (is.null(prior.sigma))
    prior.sigma <- diag(100, dim, dim)
  if (is.null(nchains))
    nchains <- 2 * dim
  S.prop <- diag(V.proposal, dim, dim)
  Theta <- array(NA, c(main.iters, dim, nchains))
  if (is.null(startVals)) {
    suppressMessages(mple <- ergm(formula, estimate = "MPLE",
                                  verbose = FALSE,
                                  offset.coef = offset.coef)$coef)
    theta <- matrix(mple + runif(dim * nchains, min = -0.1,
                                 max = 0.1), dim, nchains)
  } else {
    theta <- matrix(startVals + runif(dim * nchains, min = -0.1,
                                      max = 0.1), dim, nchains)
  }
  theta[model$etamap$offsettheta,] <- offset.coef
  acc.counts <- rep(0L, nchains)
  theta1 <- rep(NA, dim)
  tot.iters <- burn.in + main.iters
  impNet <- y
  f <- as.character(formula)
  currentFormula <- formula(paste("impNet", f[3:length(f)],
                                  sep = " ~ "))
  clock.start <- Sys.time()
  message(" > MCMC start")
  for (k in 1:tot.iters) {
    for (h in 1:nchains) {
      theta1 <- theta[, h] + gamma * apply(theta[, sample(seq(1,
                  nchains)[-h], 2)], 1, diff) + rmvnorm(1, sigma = S.prop)[1,]

      theta1[model$etamap$offsettheta] <- offset.coef

      delta <- ergm_MCMC_slave(Clist = Clist, proposal = proposal,
                               eta = theta1, control = control,
                               verbose = FALSE)$s
      pr <- dmvnorm(rbind(theta1, theta[, h]), mean = prior.mean,
                    sigma = prior.sigma, log = TRUE)
      beta <- (theta[, h] - theta1) %*% t(delta) + pr[1] -
        pr[2]
      if (beta >= log(runif(1))) {
        theta[, h] <- theta1
        if (k > burn.in) {
          acc.counts[h] <- acc.counts[h] + 1
        }
        impNet <- simulate(currentFormula, coef = theta1,
                           output = "network", basis = y, constraints = ~fixallbut(missingTies),
                           nsim = 1, control = control.simulate(MCMC.burnin = missingUpdate))
        y2 <- ergm.getnetwork(currentFormula)
        model2 <- ergm_model(currentFormula, y2)
        Clist <- ergm.Cprepare(y2, model2)
      }
    }
    if (k > burn.in)
      Theta[k - burn.in, , ] <- theta
    if (!is.null(nImp)) {
      if ((k - burn.in) == impIter * thinImp) {
        impNets[[impIter]] <- impNet
        impIter <- impIter + 1
      }
    }
  }
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)
  FF <- mcmc(apply(Theta, 2, cbind))
  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL
  ess <- round(effectiveSize(FF), 0)
  names(ess) <- model$coef.names
  out = list(Time = runtime, formula = formula, specs = model$coef.names,
             dim = dim, Theta = FF, AR = AR, ess = ess, impNets = impNets)
  class(out) <- "bergm"
  return(out)
}
