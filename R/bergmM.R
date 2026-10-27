#' Parameter estimation for Bayesian ERGMs under missing data
#'
#' Function to fit Bayesian exponential random graphs models under missing data
#' using the approximate exchange algorithm with Bayesian data augmentation
#' (Koskinen, Robins and Pattison, 2010).
#'
#' Every chain of the population carries its own imputed network, and
#' missing tie variables are updated by a Gibbs-type step from the ERGM
#' restricted to the missing dyads given the chain's current parameter.
#' This is a valid data-augmentation sampler for each chain (in Bergm 5 a
#' single imputed network was shared by all chains). Imputation steps act
#' directly on \pkg{ergm}'s compiled network states, which makes imputing at
#' every iteration cheap; network objects are only built for the returned
#' imputations.
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
#' @param prior.sigma square matrix (or scalar/vector diagonal); variance/covariance matrix for the multivariate Normal prior.
#' By default set to a diagonal matrix with every diagonal entry equal to 100.
#'
#' @param nchains count; number of chains of the population MCMC.
#' By default twice the number of free model parameters (at least 4).
#'
#' @param gamma scalar; parallel adaptive direction sampling move factor
#' (starting value when \code{adapt = TRUE}).
#'
#' @param V.proposal count; diagonal entry for the multivariate Normal proposal.
#' By default set to 0.0025.
#'
#' @param seed count; random number seed for the Bergm estimation.
#'
#' @param startVals vector; optional starting values for the parameter estimation.
#'
#' @param offset.coef vector; a vector of coefficients for the offset terms.
#'
#' @param constraints formula; constraints on the support of the network
#' distribution, as in \code{\link[ergm]{ergm}}.
#'
#' @param thin count; the thinning interval between consecutive observations.
#'
#' @param saveEveryX count; if not \code{NULL}, the partial posterior sample and imputations
#' are saved to \code{saveEveryXName} every \code{saveEveryX} iterations.
#'
#' @param cut.reject logical; if \code{TRUE}, repeated (rejected) draws are removed from
#' the pooled sample (backward compatibility only; default \code{FALSE}).
#'
#' @param saveEveryXName character; the name of the partial estimation object.
#'
#' @param imputeAllItr logical; default \code{TRUE}: missing data are imputed for every
#' chain at every iteration (exact data augmentation). If \code{FALSE}, a chain's
#' imputation is only refreshed after an accepted move, as in Bergm 5.
#'
#' @param imputeLast logical; kept for backward compatibility. Imputations always use the
#' chain's current parameter value, which is what makes the sampler valid.
#'
#' @param nImp count; number of imputed networks to be returned. If \code{NULL}, no imputed
#' network is returned. Imputations are taken at evenly spaced post burn-in
#' iterations, cycling through the chains.
#'
#' @param missingUpdate count; number of tie updates in each imputation step.
#' By default equal to the number of missing ties.
#'
#' @param imputeData data.frame; a data.frame with all attribute variables that should be
#' imputed and additional attributes that should be used for the imputation.
#' All non-numeric variables need to be specified as factors.
#' Names of vertex attributes and variable names need to be identical.
#'
#' @param attributeNames character vector; the names of all variables that need to be imputed.
#'
#' @param miceIt count; number of iterations in the MICE imputation. Default is 5.
#'
#' @param onlyKeepImputation logical; should only imputations be returned, and no bergm estimate.
#'
#' @param adapt logical; tune \code{gamma} during burn-in (see \code{\link{bergm}}).
#'
#' @param target.accept scalar; target acceptance rate for \code{adapt}.
#'
#' @param log.prior function; optional user-defined log prior (see \code{\link{bergm}}).
#'
#' @param prior.decay named vector or list with elements \code{meanlog} and
#' \code{sdlog}; log-normal prior for the decay parameters of curved terms
#' (e.g. \code{gwesp.decay}, \code{gwdegree.decay}, \code{altkstar.lambda}),
#' recycled over all of them or given one value per decay parameter. Default:
#' median 0.5, \code{sdlog} 0.75. Decay parameters are sampled on the log
#' scale, so they are always positive; \code{prior.mean} and
#' \code{prior.sigma} apply to the other parameters and may be given either
#' for all parameters (decay entries are ignored) or for the non-decay ones only.
#'
#' @param verbose logical; print progress information.
#'
#' @param ... additional arguments passed to \code{\link[ergm]{control.simulate.formula}}.
#'
#' @references
#' Caimo, A. and Friel, N. (2011), "Bayesian Inference for Exponential Random Graph Models,"
#' Social Networks, 33(1), 41-55. \url{https://arxiv.org/abs/1007.5192}
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
                   constraints = NULL,
                   thin = 1,
                   saveEveryX = NULL,
                   cut.reject = FALSE,
                   saveEveryXName = "partialBergmEstimate.rda",
                   imputeAllItr = TRUE,
                   imputeLast = TRUE,
                   nImp = NULL,
                   missingUpdate = NULL,
                   imputeData = NULL,
                   attributeNames = NULL,
                   miceIt = 5,
                   onlyKeepImputation = FALSE,
                   adapt = TRUE,
                   target.accept = 0.25,
                   log.prior = NULL,
                   prior.decay = c(meanlog = log(0.5), sdlog = 0.75),
                   verbose = TRUE,
                   ...) {
  clock.start <- Sys.time()
  this.call <- match.call()
  .check_iters(main.iters, thin)
  .with_seed(seed)
  control.args <- list(...)
  if (!isTRUE(imputeLast))
    message("'imputeLast = FALSE' is ignored: imputations use each chain's current parameter.")

  y <- ergm.getnetwork(formula)
  am <- .adjacency(y)
  has_missing <- anyNA(am)
  has_attr <- !is.null(imputeData)
  if (has_attr && is.null(attributeNames))
    stop("'attributeNames' must be given together with 'imputeData'.", call. = FALSE)
  if (!has_missing && !has_attr)
    message("Network has no missing data and no attribute data to impute were given; ",
            "no imputation will be provided.")

  imputeData2 <- NULL
  if (has_attr) {
    imputeData2 <- .imputeAttributes(y, attributeNames, imputeData, miceIt)
    y <- .set_attributes(y, imputeData2, attributeNames)
  }

  missing_el <- NULL
  if (has_missing) {
    idx <- which(is.na(am), arr.ind = TRUE)
    if (!is.directed(y)) idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
    missing_el <- as.edgelist(idx, n = network.size(y), directed = is.directed(y),
                              bipartite = if (is.bipartite(y)) y %n% "bipartite" else FALSE)
  }
  if (is.null(missingUpdate)) missingUpdate <- max(1L, sum(is.na(am)))

  build <- function(net) .bergmM_states(formula, net, constraints, missing_el, control.args)
  base <- build(y)
  setup <- base$setup
  p <- setup$p
  .check_offsets(setup, offset.coef)
  fml <- setup$formula

  mp <- .mple(fml, setup, offset.coef, constraints)
  mple <- mp$coef
  fx <- .fix_infinite(mple, setup$offset, offset.coef, verbose)
  offset <- fx$offset
  offset.coef <- fx$offset.coef
  free <- !offset
  pf <- sum(free)
  dec <- .decay_index(setup, offset)
  pm <- .param_map(setup, dec)
  prior <- .make_prior(prior.mean, prior.sigma, p, free, setup$coef.names, log.prior,
                       dec = dec, prior.decay = prior.decay)
  nchains <- as.integer(nchains %||% max(2L * pf, 4L))
  if (nchains < 3L) stop("'nchains' must be at least 3.", call. = FALSE)
  theta0 <- .initial_population(p, nchains, free, offset.coef, mple, mp$vcov,
                                startVals, setup$coef.names, offset,
                                .valid_logprior(prior$logdens, pm$eta),
                                dec = dec, decay = prior$decay)
  U_noise <- .chol_pd(.as_cov(V.proposal, p, "V.proposal")[free, free, drop = FALSE])
  eta_fn <- pm$eta

  chains <- rep(list(base[c("main", "imp", "attr")]), nchains)
  for (h in seq_len(nchains)) chains[[h]]$attr <- imputeData2

  impute_one <- function(h, theta_h, accepted) {
    ch <- chains[[h]]
    if (!(accepted || imputeAllItr)) return(as.numeric(ch$main$stats))
    if (has_missing) {
      z <- .simulate_state(ergm_state_send(ch$imp), eta_fn(theta_h), burnin = missingUpdate)
      if (!is.null(z)) {
        ch$imp <- update(ch$imp, state = z$state)
        ch$main <- update(ch$main, state = z$state[c("el", "stats", "ext.state", "ext.flag")])
      }
    }
    if (has_attr) {
      net <- as.network(ch$main)
      ch$attr <- .imputeAttributes(net, attributeNames, imputeData, miceIt)
      net <- .set_attributes(net, ch$attr, attributeNames)
      nb <- build(net)
      ch$main <- nb$main
      ch$imp <- nb$imp
    }
    chains[[h]] <<- ch
    as.numeric(ch$main$stats)
  }

  sim_one <- function(h, eta) {
    z <- .simulate_state(ergm_state_send(chains[[h]]$main), eta, aux.iters)
    if (is.null(z)) NULL else z$s[1L, ]
  }

  impNets <- impAttr <- NULL
  nImp <- if (is.null(nImp)) 0L else max(0L, min(as.integer(nImp), main.iters))
  if (nImp > 0L) {
    impNets <- if (has_missing) vector("list", nImp)
    impAttr <- if (has_attr) vector("list", nImp)
    imp_at <- burn.in + as.integer(main.iters / nImp) * seq_len(nImp)
    imp_count <- 0L
  }
  do_impute <- has_missing || has_attr
  wrapped_impute <- if (do_impute) function(h, theta_h, accepted) {
    out <- impute_one(h, theta_h, accepted)
    if (nImp > 0L && h == nchains) collect()
    out
  }
  iter <- 0L
  collect <- function() {
    iter <<- iter + 1L
    if (iter %in% imp_at) {
      imp_count <<- imp_count + 1L
      h <- (imp_count - 1L) %% nchains + 1L
      if (has_missing) impNets[[imp_count]] <<- as.network(chains[[h]]$main)
      if (has_attr) impAttr[[imp_count]] <<- chains[[h]]$attr
    }
  }

  checkpoint <- NULL
  if (!is.null(saveEveryX))
    checkpoint <- list(at = .checkpoint_at(saveEveryX, burn.in + main.iters),
                       file = saveEveryXName,
                       extra = function() list(impNets = impNets, impAttr = impAttr))

  if (verbose) message(" > MCMC start")
  res <- .exchange_engine(theta = theta0, free = free, eta_fn = eta_fn, out_fn = pm$to_nat,
                          logprior = prior$logdens,
                          sy = matrix(setup$sy, setup$q, nchains),
                          sim_one = sim_one, burn.in = burn.in, main.iters = main.iters,
                          thin = thin, gamma = gamma, U_noise = U_noise, adapt = adapt,
                          target.accept = target.accept,
                          impute = wrapped_impute, checkpoint = checkpoint,
                          blocked = FALSE, verbose = verbose)
  runtime <- difftime(Sys.time(), clock.start)

  if (onlyKeepImputation) return(list(impNets = impNets, impAttr = impAttr))
  .bergm_output(res, setup, offset, prior, main.iters, thin, cut.reject, fml, runtime,
                extra = list(mple = mple, call = this.call, impNets = impNets,
                             impAttr = impAttr,
                             settings = list(burn.in = burn.in, main.iters = main.iters,
                                             aux.iters = aux.iters, thin = thin,
                                             missingUpdate = missingUpdate,
                                             imputeAllItr = imputeAllItr)))
}

.bergmM_states <- function(formula, net, constraints, missing_el, control.args) {
  setup <- .ergm_setup(formula, constraints, control.args, network = net)
  imp <- NULL
  if (!is.null(missing_el)) {
    user <- if (is.null(constraints)) NULL else
      paste(deparse(constraints[[length(constraints)]]), collapse = " ")
    rhs <- "fixallbut(.bergm_missing)"
    if (!is.null(user) && !identical(user, ".")) rhs <- paste(rhs, "+", user)
    impc <- stats::as.formula(paste("~", rhs))
    environment(impc) <- list2env(list(.bergm_missing = missing_el),
                                  parent = environment(constraints %||% formula))
    ctrl <- do.call(control.simulate.formula, c(list(MCMC.burnin = 1L, MCMC.interval = 1L),
                                                .term_opts(control.args, network.size(net))))
    imp <- suppressMessages(
      simulate(setup$formula, coef = rep(0, setup$p), nsim = 1, constraints = impc,
               control = ctrl, return.args = "ergm_state"))$object
  }
  list(setup = setup, main = setup$state, imp = imp)
}

.set_attributes <- function(net, data, attributeNames) {
  for (i in attributeNames) {
    v <- data[, i]
    if (is.factor(v)) v <- as.character(v)
    set.vertex.attribute(net, i, v)
  }
  net
}

.imputeAttributes <- function(y, attributeNames, imputeData, miceIt = 5) {
  imputeData2 <- imputeData
  y2n <- .adjacency(y)
  mode_among <- function(v, alters) {
    tb <- table(v[alters])
    if (length(tb) == 0L || all(tb == 0)) NA_character_ else names(which.max(tb))
  }
  n <- nrow(y2n)
  impMat <- list()
  directed <- is.directed(y)
  if (directed) {
    imputeData2$indegreeImp <- colSums(y2n, na.rm = TRUE)
    imputeData2$outdegreeImp <- rowSums(y2n)
  } else {
    imputeData2$degreeImp <- rowSums(y2n)
  }
  for (a in attributeNames) {
    v <- imputeData[, a]
    if (is.numeric(v)) {
      if (directed) {
        inA <- colSums(y2n * v, na.rm = TRUE) / colSums(y2n, na.rm = TRUE)
        outA <- rowSums(sweep(y2n, 2, v, "*"), na.rm = TRUE) / rowSums(y2n, na.rm = TRUE)
        inA[is.nan(inA)] <- NA
        outA[is.nan(outA)] <- NA
        impMat[[length(impMat) + 1L]] <- inA
        impMat[[length(impMat) + 1L]] <- outA
      } else {
        avg <- rowSums(sweep(y2n, 2, v, "*"), na.rm = TRUE) / rowSums(y2n, na.rm = TRUE)
        avg[is.nan(avg)] <- NA
        impMat[[length(impMat) + 1L]] <- avg
      }
    }
    inMax <- vapply(seq_len(n), function(j) mode_among(v, which(y2n[, j] %in% 1)), "")
    impMat[[length(impMat) + 1L]] <- as.factor(inMax)
    if (directed) {
      outMax <- vapply(seq_len(n), function(j) mode_among(v, which(y2n[j, ] %in% 1)), "")
      impMat[[length(impMat) + 1L]] <- as.factor(outMax)
    }
  }
  if (length(impMat)) {
    impMat <- as.data.frame(impMat)
    names(impMat) <- paste0("impVarLongNameNoOneWillUse", seq_along(impMat))
    imputeData2 <- cbind(imputeData2, impMat)
  }
  imputeData2 <- mice::complete(mice::mice(imputeData2, m = 1, printFlag = FALSE,
                                           maxit = miceIt, remove_collinear = FALSE))
  imputeData2[, names(imputeData2) %in% names(imputeData), drop = FALSE]
}
