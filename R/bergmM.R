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
#' By default set to a vector of 0's. Note that several ergm.terms add more than one parameter to the model. You need to adjust your priors accordingly
#'
#' @param prior.sigma square matrix; variance/covariance matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100. Note that several ergm.terms add more than one parameter to the model. You need to adjust your priors accordingly
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
#' a vector of coefficients for the offset terms.
#'
#' @param constraints formula;
#' A formula specifying one or more constraints on the support of the distribution of the networks being modeled, using syntax similar to the formula argument, on the right-hand side. Multiple constraints may be given, separated by “+” and “-” operators. (See \code{\link[ergm]{ergm}} constraints for the explanation of their semantics.) Together with the model terms in the formula and the reference measure, the constraints define the distribution of networks being modeled.
#'
#' @param thin count;
#' The thinning interval between consecutive observations.
#'
#' @param cut.reject logical;
#' default TRUE. By default bergm()/bergmM() will save the last accepted theta to the posterior if the new proposal is rejected.
#' This artificially will increase the auto-correlation of consecutive parameters.
#'
#' @param saveEveryX count; If not NULL, the posterior and obtained imputation (only if nImp > 0) will be saved to your working directory at every X iterations.
#' Note that this slows down estimation and thus X should not be set too low. By default, the saved data will be in 'partialBergmEstimate.rda' and will be overwritten every X iterations and by other calls of bergm() or bergmM().
#'
#' @param saveEveryXName character; the Name of the partial estimation object.
#' If you run multiple bergm()/bergmM() calls in the same working directory you should change this name so that the calls do not override each other.
#'
#' @param imputeAllItr logical;
#' default FALSE. If TRUE, missing network and attribute data are imputed after every iteration.
#' This leads to much (!!!) longer estimation times, but potentially overall better estimation.
#' It is recommended to initially estimate with imputeAllItr = FALSE to get reasonable starting values and reduce overall estimation time.
#'
#' @param imputeLast logical;
#' default TRUE and in line with Koskinen et al. 2010.
#' If FALSE, network imputations will be performed with a random draw from the so far accepted parameter values.
#'
#' @param nImp count;
#' number of imputed networks to be returned. If null, no imputed network will be returned.
#'
#' @param missingUpdate count;
#' number of tie updates in each imputation step.
#' By default equal to the number of missing ties.
#' Smaller numbers increase speed. Larger numbers lead to better sampling.
#'
#' @param imputeData data.frame;
#' a data.frame with all attribute variables that should be imputed and additional attributes that should be used for the imputation.
#' All non-numeric variables need to be specified as.factors.
#' Names of vertex.attributes and variable names need to be identical.
#'
#' @param attributeNames character vector,
#' a vector with the names of all variables that need to be imputed.
#' These names need to be identical with the vertex.attributes and must be part of the names of the \code{imputeData} data.frame.
#'
#' @param miceIt count,
#' number of iterations in the MICE imputation. Default is 5.
#'
#' @param onlyKeepImputation logical,
#' Should only imputations be returned, and no bergm estimate (only recommended after you made sure that the model estimates properly).
#'
#' @param ... additional arguments, to be passed to lower-level functions.
#'
#' @references
#' Caimo, A. and Friel, N. (2011), "Bayesian Inference for Exponential Random Graph Models,"
#' Social Networks, 33(1), 41-55. \url{https://arxiv.org/abs/1007.5192}
#'
#' Caimo, A. and Friel, N. (2014), "Bergm: Bayesian Exponential Random Graphs in R,"
#' Journal of Statistical Software, 61(2), 1-25. \url{https://www.jstatsoft.org/v61/i02}
#'
#' Koskinen, J.H., Robins, G.L., and Pattison, P.E. (2010), "Analysing exponential
#' random graph (p-star) models with missing data using Bayesian data augmentation,"
#' Statistical Methodology 7(3), 366-384.
#'
#' Krause, R.W., Huisman, M., Steglich, C., and Snijders, T.A. (2020), "Missing data in
#' cross-sectional networks-An extensive comparison of missing data treatment methods",
#' Social Networks 62: 99-112.
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network
#' data(florentine)
#'
#' # Create missing data
#' set.seed(14021994)
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
bergmM <- function (formula, burn.in = 100, main.iters = 1000, aux.iters = 1000, 
          prior.mean = NULL, prior.sigma = NULL, nchains = NULL, gamma = 0.5, 
          V.proposal = 0.0025, seed = NULL, startVals = NULL, offset.coef = NULL, 
          constraints = NULL, thin = 1, saveEveryX = NULL, saveEveryXName = "partialBergmEstimate.rda", 
          imputeAllItr = FALSE, imputeLast = TRUE, nImp = NULL, missingUpdate = NULL, 
          imputeData = NULL, attributeNames = NULL, miceIt = 5, onlyKeepImputation = FALSE, 
          ...) 
{
  if (is.null(seed)) {
    set.seed(sample(1:999, 1))
  }
  else {
    set.seed(seed)
  }
  y <- ergm.getnetwork(formula)
  imputeAttributes <- function(y, attributeNames, imputeData) {
    imputeData2 <- imputeData
    y2n <- as.matrix.network(y)
    if (y$gal$directed) {
      imputeData2$indegreeImp <- colSums(y2n, na.rm = TRUE)
      imputeData2$outdegreeImp <- rowSums(y2n)
      impMat <- as.data.frame(matrix(NA, nrow = nrow(y2n), 
                                     ncol = 0))
      for (i in 1:length(attributeNames)) {
        if (is.numeric(imputeData[, attributeNames[i]])) {
          avgInAlt <- rowSums(sweep(t(y2n), MARGIN = 2, 
                                    imputeData[, attributeNames[i]], "*"), na.rm = TRUE)/rowSums(t(y2n), 
                                                                                                 na.rm = TRUE)
          avgInAlt[is.nan(avgInAlt)] <- NA
          impMat <- cbind(impMat, avgInAlt)
          avgOutAlt <- rowSums(sweep(y2n, MARGIN = 2, 
                                     imputeData[, attributeNames[i]], "*"), na.rm = TRUE)/rowSums(y2n, 
                                                                                                  na.rm = TRUE)
          avgOutAlt[is.nan(avgOutAlt)] <- NA
          impMat <- cbind(impMat, avgOutAlt)
        }
        inMax <- c()
        outMax <- c()
        for (j in 1:nrow(y2n)) {
          inMax <- c(inMax, ifelse(is.null(names(which.max(table(imputeData[, 
                                                                            attributeNames[i]][as.logical(y2n[, j])])))), 
                                   yes = NA, no = names(which.max(table(imputeData[, 
                                                                                   attributeNames[i]][as.logical(y2n[, j])])))))
          outMax <- c(outMax, ifelse(is.null(names(which.max(table(imputeData[, 
                                                                              attributeNames[i]][as.logical(y2n[j, ])])))), 
                                     yes = NA, no = names(which.max(table(imputeData[, 
                                                                                     attributeNames[i]][as.logical(y2n[j, ])])))))
        }
        impMat <- cbind(impMat, as.factor(inMax))
        impMat <- cbind(impMat, as.factor(outMax))
      }
      names(impMat) <- c(paste("impVarLongNameNoOneWillUse", 
                               1:ncol(impMat), sep = ""))
      imputeData2 <- cbind(imputeData2, impMat)
    }
    else {
      imputeData2$degreeImp <- rowSums(y2n)
      impMat <- as.data.frame(matrix(NA, nrow = nrow(y2n), 
                                     ncol = 0))
      for (i in 1:length(attributeNames)) {
        if (is.numeric(imputeData[, attributeNames[i]])) {
          avgAlt <- rowSums(sweep(y2n, MARGIN = 2, imputeData[, 
                                                              attributeNames[i]], "*"), na.rm = TRUE)/rowSums(y2n, 
                                                                                                              na.rm = TRUE)
          avgAlt[is.nan(avgAlt)] <- NA
          impMat <- cbind(impMat, avgAlt)
        }
        inMax <- c()
        for (j in 1:nrow(y2n)) {
          inMax <- c(inMax, ifelse(is.null(names(which.max(table(imputeData[, 
                                                                            attributeNames[i]][as.logical(y2n[, j])])))), 
                                   yes = NA, no = names(which.max(table(imputeData[, 
                                                                                   attributeNames[i]][as.logical(y2n[, j])])))))
        }
        impMat <- cbind(impMat, as.factor(inMax))
      }
      names(impMat) <- c(paste("impVarLongNameNoOneWillUse", 
                               1:ncol(impMat), sep = ""))
      imputeData2 <- cbind(imputeData2, impMat)
    }
    imputeData2 <- complete(mice(imputeData2, m = 1, printFlag = FALSE, 
                                 maxit = miceIt, remove_collinear = FALSE))
    imputeData2 <- imputeData2[, names(imputeData2) %in% 
                                 names(imputeData)]
    return(imputeData2)
  }
  if (!is.null(imputeData)) {
    imputeData2 <- imputeAttributes(y = y, attributeNames = attributeNames, 
                                    imputeData = imputeData)
    for (i in attributeNames) {
      if (is.factor(imputeData2[, i])) {
        imputeData2[, i] <- as.character(imputeData2[, 
                                                     i])
      }
      set.vertex.attribute(y, i, imputeData2[, i])
    }
  }
  model <- ergm_model(formula, y)
  specs <- unlist(sapply(model$terms, "[", "coef.names"), use.names = FALSE)
  if (!is.null(imputeData)) {
    formula <- as.formula(paste("y", "~", as.character(formula)[3]))
  }
  sy <- summary(formula)
  dim <- length(sy)
  if (dim == 1) {
    stop("Model dimension must be greater than 1")
  }

  if (!any(is.na(as.matrix.network(y))) && is.null(imputeData)) {
    print("Network has no missing data. \n\n          No attribute data to impute was given. \n\n          No imputation will be provided.")
  }
  impNets <- NULL
  if (!is.null(nImp)) {
    nImp <- max(0, min(nImp, main.iters))
    thinImp <- as.integer(main.iters/nImp)
    impIter <- 1
    impNets <- vector("list", nImp)
    if (!is.null(imputeData)) {
      impAttr <- vector("list", nImp)
    }
  }
  missingTies <- matrix(0, y$gal$n, y$gal$n)
  missingTies[is.na(as.matrix.network(y))] <- 1
  missingTies <- as.edgelist(as.network(missingTies), n = y$gal$n)
  if (is.null(missingUpdate)) {
    missingUpdate <- sum(is.na(as.matrix.network(y)))
  }
  impNet <- y
  f <- as.character(formula)
  currentFormula <- formula(paste("impNet", f[3:length(f)], 
                                  sep = " ~ "))
  if (is.null(constraints)) {
    constraints <- ~.
    impConstraints <- as.formula("~ fixallbut(missingTies)")
  } else {
    impConstraints <- as.formula(str_c("~ fixallbut(missingTies)", 
                                       str_split(constraints, pattern = "~")[[2]][1], sep = " + "))
  }
  y0 <- simulate(currentFormula, coef = rep(0, dim), nsim = 1, 
                 control = control.simulate(MCMC.burnin = 1, MCMC.interval = 1), 
                 return.args = "ergm_state", constraints = constraints)$object
  control <- control.ergm(MCMC.burnin = aux.iters, MCMC.interval = 1, 
                          MCMC.samplesize = 1, ...)
  if (!is.null(control$init)) {
    if (length(control$init) != length(model$etamap$offsettheta)) {
      stop("Invalid starting parameter vector control$init: \n\n            wrong number of parameters. \n\n            If you are passing output from another ergm run as control$init, \n\n            in a model with curved terms, see help(enformulate.curved).")
    }
  } else {
    control$init <- rep(NA, length(model$etamap$offsettheta))
  }
  if (!is.null(offset.coef)) {
    if (length(control$init[model$etamap$offsettheta]) != 
        length(offset.coef)) {
      stop("Invalid offset parameter vector offset.coef: \n\n           wrong number of parameters: expected ", 
           length(control$init[model$etamap$offsettheta]), 
           " got ", length(offset.coef), ".")
    }
    control$init[model$etamap$offsettheta] <- offset.coef
  }
  if (any(is.na(control$init) & model$etamap$offsettheta)) {
    stop("The model contains offset terms whose parameter values have not been specified:", 
         paste.and(specs[is.na(control$init) | model$offsettheta]), 
         ".", sep = "")
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
  suppressMessages(mple <- stats::coef(ergm(formula, estimate = "MPLE", 
                                            verbose = FALSE, offset.coef = offset.coef)))
  if (any(is.infinite(mple))) {
    offset.coeffs <- rep(NA, length(model$etamap$offsettheta))
    if (!is.null(offset.coef)) {
      offset.coeffs[model$etamap$offsettheta] <- offset.coef
    }
    if (any(mple == -Inf)) {
      cat("Observed statistic(s) for\n", format(names(which(mple == 
                                                              -Inf))), "\nare at their smallest attainable values.", 
          "\n Their coefficients will be fixed at -Inf.", 
          "\n No such ties will be allowed to exist.\n\n")
      offset.coeffs[!model$etamap$offsettheta][which(mple == 
                                                       -Inf)] <- -Inf
    }
    if (any(mple == Inf)) {
      cat("Observed statistic(s) for\n", format(names(which(mple == 
                                                              Inf))), "\nare at their largest attainable values.", 
          "\n Their coefficients will be fixed at Inf.", 
          "\n All such ties will be forced to exist.\n\n")
      offset.coeffs[!model$etamap$offsettheta][which(mple == 
                                                       Inf)] <- Inf
    }
    model$etamap$offsettheta[!model$etamap$offsettheta][which(is.infinite(mple))] <- TRUE
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
  acc.counts <- rep(0L, nchains)
  theta1 <- rep(NA, dim)
  tot.iters <- burn.in + main.iters
  lastPar <- NULL
  

  clock.start <- Sys.time()
  message(" > MCMC start")
  for (k in 1:tot.iters) {
    for (h in 1:nchains) {
      accepted <- FALSE
      theta1 <- theta[, h] + 
        gamma * apply(theta[, sample(seq(1, nchains)[-h],2)], 1, diff) + 
        rmvnorm(1, sigma = S.prop)[1, ]
      theta1[model$etamap$offsettheta] <- offset.coef
      delta <- ergm_MCMC_sample(y0, 
                                theta = theta1,
                                control = control)$stats[[1]][1, ] - sy
      pr <- dmvnorm(rbind(theta1[!model$etamap$offsettheta], 
                          theta[!model$etamap$offsettheta, h]), 
                    mean = prior.mean[!model$etamap$offsettheta], 
                    sigma = prior.sigma[!model$etamap$offsettheta,
                                        !model$etamap$offsettheta],
                    log = TRUE)
      beta <- (theta[!model$etamap$offsettheta, h] - 
                 theta1[!model$etamap$offsettheta]) %*% 
        delta[!model$etamap$offsettheta] + pr[1] - pr[2]
      if (beta >= log(runif(1))) {
        theta[, h] <- theta1
        if (k > burn.in) {
          Theta[k - burn.in, , ] <- theta
        }
        accepted <- TRUE
        lastPar <- theta1
        if (k > burn.in) {
          acc.counts[h] <- acc.counts[h] + 1
        }
      }
      if (accepted || (imputeAllItr && !is.null(lastPar))) {
        if (any(is.na(as.matrix.network(y)))) {
          if (!all(is.na(Theta)) && !imputeLast) {
            impPars <- Theta[sample(1:(k - burn.in), 
                                    1), , sample(1:nchains, 1)]
            counter <- 0
            while (any(is.na(impPars))) {
              impPars <- Theta[sample(1:(k - burn.in), 
                                      1), , sample(1:nchains, 1)]
              counter <- counter + 1
              if (counter == nchains * 2) {
                impPars <- lastPar
              }
            }
          } else {
            impPars <- lastPar
          }
          impNet <- simulate(currentFormula, coef = impPars, 
                             output = "network", 
                             basis = impNet, 
                             constraints = impConstraints, 
                             nsim = 1, 
                             control = control.simulate(MCMC.burnin = missingUpdate))
        }
        if (!is.null(imputeData)) {
          imputeData2 <- imputeAttributes(y = impNet, 
                                          attributeNames = attributeNames, 
                                          imputeData = imputeData)
          for (i in attributeNames) {
            if (is.factor(imputeData2[, i])) {
              imputeData2[, i] <- as.character(imputeData2[,i])
            }
            set.vertex.attribute(impNet, i, imputeData2[,i])
          }
        }
        y0 <- simulate(currentFormula,
                       coef = rep(0, dim),
                       nsim = 1, 
                       control = control.simulate(MCMC.burnin = 1,
                                                  MCMC.interval = 1), 
                       return.args = "ergm_state", 
                       constraints = constraints)$object
        sy <- summary(currentFormula)
      }
    }
    if (k > burn.in) {
      Theta[k - burn.in, , ] <- theta
    }
    if (!is.null(nImp)) {
      if ((k - burn.in) == impIter * thinImp) {
        if (any(is.na(as.matrix.network(y)))) {
          impNets[[impIter]] <- impNet
        }
        if (!is.null(imputeData)) {
          impAttr[[impIter]] <- imputeData2
        }
        impIter <- impIter + 1
      }
    }
    if (!is.null(saveEveryX)) {
      if (k %in% saveXseq) {
        if (is.null(nImp)) {
          impNets <- NULL
        }
        if (is.null(imputeData)) {
          impAttr <- NULL
        }
        partialBergmEstimate <- list(Theta = Theta, impNets = impNets, 
                                     impAttr = impAttr, failedAfter = k)
        save(partialBergmEstimate, file = saveEveryXName)
      }
    }
  }
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)
  clock.end <- Sys.time()
  runtime <- difftime(clock.end, clock.start)
  Theta <- apply(Theta, 2, cbind)
  FF <- mcmc(Theta)
  colnames(FF) <- names(mple)
  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL
 FF <- mcmc(Theta)
  if (cut.reject) {
    FF <- as.matrix(unique(as.data.frame(FF)))
  }
  if (thin > 1) {
    FF <- FF[seq(1, nrow(FF), thin), ]
  }
  colnames(FF) <- names(mple)
  ess_out <- rep(NA, dim)
  ess <- round(effectiveSize(FF[, !model$etamap$offsettheta]), 
               0)
  ess_out[!model$etamap$offsettheta] <- ess
  names(ess_out) <- names(mple)
  fixed <- model$etamap$offsettheta
  names(fixed) <- names(mple)
  class(FF) <- 'mcmc'
  AR <- round(1 - rejectionRate(FF)[1], 2)
  names(AR) <- NULL
 
  if (is.null(imputeData)) {
    impAttr <- NULL
  }
  if (onlyKeepImputation) {
    out <- list(impNets = impNets, impAttr = impAttr)
    return(out)
  }
  out <- list(Time = runtime, formula = formula, specs = specs, 
              dim = dim, Theta = mcmc(unique(as.matrix(FF))), AR = AR, 
              ess = ess, impNets = impNets, impAttr = impAttr)
  class(out) <- "bergm"
  return(out)
}
