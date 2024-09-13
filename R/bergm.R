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
#' @param offset.coef vector;
#' A vector of coefficients for the offset terms.
#'
#' @param saveEveryX count; If not NULL, the posterior and obtained imputation (only if nImp > 0) will be saved to your working directory at every X iterations.
#' Note that this slows down estimation and thus X should not be set too low. By default, the saved data will be in 'partialBergmEstimate.rda' and will be overwritten every X iterations and by other calls of bergm() or bergmM().
#'
#' @param saveEveryXName character; the Name of the partial estimation object.
#' If you run multiple bergm()/bergmM() calls in the same working directory you should change this name so that the calls do not override each other.
#'
#' @param ... additional arguments, to be passed to lower-level functions.
#'
#' @references
#' Caimo, A. and Friel, N. (2011), "Bayesian Inference for Exponential Random Graph Models,"
#' Social Networks, 33(1), 41-55. \url{https://arxiv.org/abs/1007.5192}
#'
#' Caimo, A. and Friel, N. (2014), "Bergm: Bayesian Exponential Random Graphs in R,"
#' Journal of Statistical Software, 61(2), 1-25. \url{https://www.jstatsoft.org/article/view/v061i02}
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data(florentine)
#'
#' # Posterior parameter estimation:
#' p.flo <- bergm(flomarriage ~ edges + kstar(2),
#'                burn.in    = 50,
#'                aux.iters  = 500,
#'                main.iters = 3000,
#'                gamma      = 1.2)
#'
#' # Posterior summaries:
#' summary(p.flo)
#' }
#' @export
#'
bergm <- function(formula,
                  prior.mean = NULL,
                  prior.sigma = NULL,
                  burn.in = 100,
                  main.iters = 1000,
                  aux.iters = 1000,
                  nchains = NULL,
                  gamma = 0.5,
                  V.proposal = 0.0025,
                  startVals = NULL,
                  offset.coef = NULL,
                  constraints = NULL,
                  thin = 1,
                  cut.reject = FALSE,
                  saveEveryX = NULL,
                  saveEveryXName = 'partialBergmEstimate.rda',
                  ...) {
  library(statnet)
  library(coda)
  library(mvtnorm)
  y <- ergm.getnetwork(formula)
  model <- ergm_model(formula, y)
  specs <- unlist(sapply(model$terms, '[', 'coef.names'), use.names = FALSE)
  sy <- summary(formula)
  dim <- length(sy)
  if (dim == 1)
    stop("Model dimension must be greater than 1.")
  if (any(is.na(as.matrix.network(y))))
    print("Network has missing data. Use bergmM() instead.")
  if (!is.null(offset.coef)) {
    if (any(offset.coef %in% c(NaN, NA))) {
      stop("NaN, NA are not allowed in offset.coef.")
    }
  }

  if (is.null(constraints)) {
    constraints <- ~.
  }

  y0 <- simulate(formula,
                 coef = rep(0, dim),
                 nsim = 1,
                 control = control.simulate(MCMC.burnin = 1,
                                            MCMC.interval = 1),
                 return.args = "ergm_state",
                 constraints = constraints)$object

  control <- control.ergm(MCMC.burnin = aux.iters,
                          MCMC.interval = 1,
                          MCMC.samplesize = 1)
  if (!is.null(control$init)) {
    if (length(control$init) != length(model$etamap$offsettheta)) {
      stop(paste("Invalid starting parameter vector control$init:",
                 "wrong number of parameters.",
                 "If you are passing output from another ergm run as control$init,",
                 "in a model with curved terms, see help(enformulate.curved)."))
    }
  } else {
    control$init <- rep(NA, length(model$etamap$offsettheta))
  }
  if (!is.null(offset.coef)) {
    if (length(control$init[model$etamap$offsettheta]) != length(offset.coef)) {
      stop("Invalid offset parameter vector offset.coef: ",
           "wrong number of parameters: expected ",
           length(control$init[model$etamap$offsettheta]),
           " got ", length(offset.coef), ".")
    }
    control$init[model$etamap$offsettheta] <- offset.coef
  }

  if (any(is.na(control$init) & model$etamap$offsettheta)) {
    stop("The model contains offset terms whose parameter values have not been
         specified:",
         paste.and(specs[is.na(control$init) | model$offsettheta]), ".",
         sep = "")
  }


  if (!is.null(saveEveryX)) {
    saveXseq <- seq(1, main.iters, saveEveryX)
  }

  if (is.null(prior.mean)) {
    prior.mean <- rep(0, dim)
  }
  if (is.null(prior.sigma)) {
    prior.sigma <- diag(100, dim, dim)
  }
  if (is.null(nchains)) {
    nchains <- 2 * dim
  }
  S.prop <- diag(V.proposal, dim, dim)
  Theta <- array(NA, c(main.iters, dim, nchains))


  suppressMessages(mple <- ergm(formula,
                                estimate = "MPLE",
                                verbose = FALSE,
                                offset.coef = offset.coef) |> stats::coef())

  if (any(is.infinite(mple))) {
    offset.coeffs <- rep(NA,length(model$etamap$offsettheta))
    if (!is.null(offset.coef)) {
      offset.coeffs[model$etamap$offsettheta] <- offset.coef
    }
    if (any(mple == -Inf)) {
      cat('Observed statistic(s) for\n',
          format(names(which(mple == -Inf))),
          '\nare at their smallest attainable values.',
          '\n Their coefficients will be fixed at -Inf.',
          '\n No such ties will be allowed to exist.\n\n')
      offset.coeffs[!model$etamap$offsettheta][which(mple == -Inf)] <- -Inf
    }
    if (any(mple == Inf)) {
      cat('Observed statistic(s) for\n',
          format(names(which(mple == Inf))),
          '\nare at their largest attainable values.',
          '\n Their coefficients will be fixed at Inf.',
          '\n All such ties will be forced to exist.\n\n')
      offset.coeffs[!model$etamap$offsettheta][which(mple == Inf)] <- Inf
    }

    model$etamap$offsettheta[!model$etamap$offsettheta][
      which(is.infinite(mple))] <- TRUE
    offset.coeffs <- na.omit(offset.coeffs)
    attributes(offset.coeffs) <- NULL
    offset.coef <- offset.coeffs
  }

  if (!is.null(startVals)) {
    theta <- matrix(startVals + runif(dim * nchains, min = -0.1,
                                      max = 0.1), dim, nchains)
  } else {
    theta <- matrix(mple + runif(dim * nchains, min = -0.1,
                                 max = 0.1), dim, nchains)
  }

  theta[model$etamap$offsettheta, ] <- offset.coef
  theta1 <- rep(NA, dim)
  tot.iters <- burn.in + main.iters
  clock.start <- Sys.time()
  prior.mean_cut <- prior.mean[!model$etamap$offsettheta]
  prior.sigma_cut <- prior.sigma[!model$etamap$offsettheta,]
  prior.sigma_cut <- prior.sigma_cut[,!model$etamap$offsettheta]
  message(" > MCMC start")
  for (k in 1:tot.iters) {
    for (h in 1:nchains) {
      theta1 <- theta[, h] +
        gamma * apply(theta[, sample(seq(1,nchains)[-h], 2)], 1, diff) +
        rmvnorm(1, sigma = S.prop)[1,]
      theta1[model$etamap$offsettheta] <- offset.coef

      delta <- ergm_MCMC_sample(y0,
                                theta = theta1,
                                control = control)$stats[[1]][1,] - sy
      delta_cut  <- delta[!model$etamap$offsettheta]
      theta1_cut <- theta1[!model$etamap$offsettheta]
      thetah_cut <- theta[!model$etamap$offsettheta,h]

      pr <- dmvnorm(rbind(theta1_cut,
                          thetah_cut),
                    mean = prior.mean_cut,
                    sigma = prior.sigma_cut,
                    log = TRUE)
      beta <- (thetah_cut - theta1_cut) %*% delta_cut + pr[1] - pr[2]
      if (beta >= log(runif(1))) {
        theta[, h] <- theta1
      }
    }
    if (k > burn.in) {
      Theta[k - burn.in, , ] <- theta
    }
    if (!is.null(saveEveryX)) {
      if (k %in% saveXseq) {
        partialBergmEstimate <- list(Theta = Theta,
                                     failedAfter = k)
        save(partialBergmEstimate, file = saveEveryXName)
      }
    }

  }
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)

  Theta <- apply(Theta, 2, cbind)
  if (cut.reject) {
    Theta <- unique(Theta)
  }

  if (thin > 1) {
    Theta <- Theta[seq(1,nrow(Theta),thin),]
  }
  FF <- mcmc(Theta)
  colnames(FF) <- names(mple)


  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL

  colnames(FF) <- names(mple)
  ess_out <- rep(NA,dim)
  ess <- round(effectiveSize(FF[,!model$etamap$offsettheta]), 0)
  ess_out[!model$etamap$offsettheta] <- ess
  names(ess_out) <- names(mple)

  fixed <- model$etamap$offsettheta
  names(fixed) <- names(mple)

  out = list(Time = runtime,
             formula = formula,
             specs = names(mple),
             dim = dim,
             Theta = FF,
             AR = AR,
             ess = ess_out,
             fixed = fixed)
  class(out) <- "bergm"
  return(out)
}
