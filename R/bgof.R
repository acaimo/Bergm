#' Bayesian goodness-of-fit diagnostics for ERGMs
#'
#' Function to calculate summaries for degree,
#' minimum geodesic distances,
#' and edge-wise shared partner distributions
#' to diagnose the Bayesian goodness-of-fit of
#' exponential random graph models.
#'
#' For every sampled posterior parameter one network is simulated from the
#' model, starting from the observed network. The goodness-of-fit statistics
#' are monitored inside \pkg{ergm}'s compiled sampler (as additional
#' zero-coefficient model terms), so no network objects are created and the
#' posterior predictive distribution of the model statistics is obtained at
#' no extra cost. Simulations can run in parallel with \code{cores > 1}.
#'
#' @param x an \code{R} object of class \code{bergm}.
#'
#' @param sample.size count; number of networks
#' to be simulated and compared to the observed network.
#'
#' @param aux.iters count; number of iterations used for network simulation.
#'
#' @param n.deg count; used to plot only the first
#' \code{n.deg}-1 degree distributions.
#' By default no restrictions on the number of degree
#' distributions is applied.
#'
#' @param n.dist count; used to plot only the first
#' \code{n.dist}-1 geodesic distances distributions.
#' By default no restrictions on the number of geodesic
#' distances distributions is applied.
#'
#' @param n.esp count; used to plot only the first
#' \code{n.esp}-1 edge-wise shared partner distributions.
#' By default no restrictions on the number of
#' edge-wise shared partner distributions is applied.
#'
#' @param n.ideg count; used to plot only the first
#' \code{n.ideg}-1 in-degree distributions.
#' By default no restrictions on the number of
#' in-degree distributions is applied.
#'
#' @param n.odeg count; used to plot only the first
#' \code{n.odeg}-1 out-degree distributions.
#' By default no restrictions on the number of
#' out-degree distributions is applied.
#'
#' @param plot logical; draw the diagnostic plots (default \code{TRUE}).
#'
#' @param cores count; number of parallel workers.
#'
#' @param seed integer; optional seed.
#'
#' @param ... additional arguments, passed to \code{\link{plot.bgof}}.
#'
#' @return An object of class \code{bgof} (returned invisibly when plotted):
#' simulated (\code{sim.*}) and observed (\code{obs.*}) proportions for each
#' goodness-of-fit statistic, plus \code{sim.model}/\code{obs.model} and
#' posterior predictive p-values (\code{pval.model}) for the model statistics.
#'
#' @references
#' Caimo, A. and Friel, N. (2011), "Bayesian Inference for Exponential Random Graph Models,"
#' Social Networks, 33(1), 41-55. \url{https://arxiv.org/abs/1007.5192}
#'
#' Caimo, A. and Friel, N. (2014), "Bergm: Bayesian Exponential Random Graphs in R,"
#' Journal of Statistical Software, 61(2), 1-25. \doi{10.18637/jss.v061.i02}
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
#'                main.iters = 1000,
#'                gamma      = 1.2)
#'
#' # Bayesian goodness-of-fit test:
#' bgof(p.flo,
#'      aux.iters   = 500,
#'      sample.size = 30,
#'      n.deg       = 10,
#'      n.dist      = 9,
#'      n.esp       = 6)
#'}
#' @export
bgof <- function(x,
                 sample.size = 100,
                 aux.iters   = 10000,
                 n.deg       = NULL,
                 n.dist      = NULL,
                 n.esp       = NULL,
                 n.ideg      = NULL,
                 n.odeg      = NULL,
                 plot        = TRUE,
                 cores       = 1,
                 seed        = NULL,
                 ...) {
  stopifnot(inherits(x, "bergm"))
  .with_seed(seed)
  Theta <- as.matrix(x$Theta)
  FF <- Theta[sample.int(nrow(Theta), sample.size, replace = sample.size > nrow(Theta)), ,
              drop = FALSE]

  pp <- .posterior_predictive(x, FF, aux.iters, gof = TRUE, cores = cores, seed = seed)
  obs <- pp$obs
  S <- pp$stats
  n <- pp$n
  grp <- pp$groups
  prop_mat <- function(cols) {
    o <- obs[cols]
    s <- S[, cols, drop = FALSE]
    po <- if (sum(o) > 0) o / sum(o) else o * 0
    rs <- rowSums(s)
    ps <- s / ifelse(rs > 0, rs, 1)
    ps[rs == 0, ] <- 1
    cbind(po, t(ps))
  }
  model_cols <- grp$model
  sim.model <- S[, model_cols, drop = FALSE]
  obs.model <- obs[model_cols]
  d <- sweep(sim.model, 2, obs.model)
  pval.model <- pmin(1, 2 * pmin(colMeans(d >= 0), colMeans(d <= 0)))

  B <- prop_mat(grp$dist)
  C <- prop_mat(grp$esp)
  out <- list(sim.dist = B[, -1, drop = FALSE], sim.esp = C[, -1, drop = FALSE],
              obs.dist = B[, 1], obs.esp = C[, 1],
              sim.model = sim.model, obs.model = obs.model, pval.model = pval.model,
              directed = pp$directed, formula = x$formula, theta = FF)
  if (!pp$directed) {
    A <- prop_mat(grp$deg)
    out <- c(list(sim.degree = A[, -1, drop = FALSE], obs.degree = A[, 1]), out)
  } else {
    A <- prop_mat(grp$ideg)
    AA <- prop_mat(grp$odeg)
    out <- c(list(sim.idegree = A[, -1, drop = FALSE], sim.odegree = AA[, -1, drop = FALSE],
                  obs.degree = A[, 1], obs.idegree = A[, 1], obs.odegree = AA[, 1]), out)
  }
  class(out) <- "bgof"
  if (plot) {
    plot.bgof(out, n.deg = n.deg, n.dist = n.dist, n.esp = n.esp,
              n.ideg = n.ideg, n.odeg = n.odeg, ...)
    return(invisible(out))
  }
  out
}

.posterior_predictive <- function(x, FF, aux.iters, gof = TRUE, cores = 1, seed = NULL,
                                  networks = FALSE) {
  formula <- x$formula
  y <- ergm.getnetwork(formula)
  n <- network.size(y)
  directed <- is.directed(y)
  if (gof && is.bipartite(y))
    stop("bgof() does not support bipartite networks; use simulate() on the fit ",
         "to compute posterior predictive statistics.", call. = FALSE)
  topts <- .term_opts(list(), n)$term.options
  model <- ergm_model(formula, y, term.options = topts)
  q_model <- nparam(model, canonical = TRUE)
  has_geodist <- exists("InitErgmTerm.geodistdist", envir = asNamespace("ergm"))
  mon <- NULL
  if (gof) {
    mon <- if (directed) {
      stats::as.formula(sprintf("~ idegree(0:%d) + odegree(0:%d) + esp(0:%d)", n - 1, n - 1, n - 2))
    } else {
      stats::as.formula(sprintf("~ degree(0:%d) + esp(0:%d)", n - 1, n - 2))
    }
    if (has_geodist) mon <- stats::update(mon, ~ . + geodistdist)
  }
  sim <- suppressMessages(simulate(formula, coef = rep(0, nparam(model, canonical = FALSE)),
                                   nsim = 1, constraints = x$constraints %||% ~.,
                                   monitor = mon, return.args = "ergm_state",
                                   control = control.simulate.formula(term.options = topts)))
  full <- sim$object
  st <- ergm_state_send(full)
  q_tot <- nparam(full$model, canonical = TRUE)
  etamap <- model$etamap
  etas <- lapply(seq_len(nrow(FF)), function(i)
    c(ergm.eta(FF[i, ], etamap), rep(0, q_tot - q_model)))

  need_states <- networks || (gof && !has_geodist)
  run <- function(eta) {
    z <- .simulate_state(st, eta, burnin = aux.iters)
    if (is.null(z)) stop(.sim_limit_msg, call. = FALSE)
    if (need_states) list(s = z$s[1L, ], state = z$state) else list(s = z$s[1L, ])
  }
  if (cores > 1 && !need_states) {
    cl <- .start_cluster(cores, seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, .worker_set, "pp_state", st)
    chunks <- .chunk(seq_along(etas), length(cl))
    res <- parallel::clusterApply(cl, lapply(chunks, function(ix) etas[ix]), .worker_sim,
                                  burnin = aux.iters, key = "pp_state")
    res <- lapply(unlist(res, recursive = FALSE), function(m) {
      if (is.null(m)) stop(.sim_limit_msg, call. = FALSE)
      list(s = m[1L, ])
    })
  } else {
    res <- lapply(etas, run)
  }
  S <- do.call(rbind, lapply(res, `[[`, "s"))
  nm <- param_names(full$model, canonical = TRUE)
  colnames(S) <- nm
  obs <- stats::setNames(as.numeric(full$stats), nm)
  nets <- NULL
  if (need_states) nets <- lapply(res, function(r) as.network(update(full, state = r$state)))

  groups <- list(model = seq_len(q_model))
  if (gof) {
    mon_idx <- seq_len(q_tot)[-seq_len(q_model)]
    pick <- function(re) mon_idx[grepl(re, nm[mon_idx])]
    if (!directed) {
      groups$deg <- pick("^degree[0-9]+$")
    } else {
      groups$ideg <- pick("^idegree[0-9]+$")
      groups$odeg <- pick("^odegree[0-9]+$")
    }
    groups$esp <- pick("^esp")
    if (has_geodist) {
      groups$dist <- pick("^geodist")
    } else {
      dd <- t(vapply(nets, ergm.geodistdist, numeric(n)))
      od <- ergm.geodistdist(y)
      k <- ncol(S)
      S <- cbind(S, dd)
      obs <- c(obs, od)
      groups$dist <- k + seq_len(n)
    }
  }
  list(stats = S, obs = obs, groups = groups, n = n, directed = directed,
       networks = if (networks) nets)
}

#' Plot Bayesian goodness-of-fit diagnostics
#'
#' @param x an object of class \code{bgof}.
#' @param n.deg,n.dist,n.esp,n.ideg,n.odeg counts; see \code{\link{bgof}}.
#' @param ... unused.
#' @export
#' @method plot bgof
plot.bgof <- function(x, n.deg = NULL, n.dist = NULL, n.esp = NULL,
                      n.ideg = NULL, n.odeg = NULL, ...) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  panel <- function(obs, sim, nshow, xlab, ylab, labels) {
    nshow <- max(1L, min(nshow, length(obs)))
    labels <- labels[seq_len(nshow)]
    keep <- seq_len(nshow)
    M <- sim[keep, , drop = FALSE]
    q5 <- apply(M, 1, stats::quantile, probs = 0.05)
    q95 <- apply(M, 1, stats::quantile, probs = 0.95)
    graphics::boxplot(as.data.frame(t(M)), xaxt = "n", xlab = xlab, ylab = ylab,
                      col = "gray90", border = "gray45", outcol = "gray45")
    graphics::axis(1, seq_len(nshow), labels)
    graphics::lines(obs[keep], lwd = 2.5, col = "black")
    graphics::lines(q5, col = "gray40", lty = 2)
    graphics::lines(q95, col = "gray40", lty = 2)
  }
  dist_panel <- function() {
    B <- cbind(x$obs.dist, x$sim.dist)
    nb <- nrow(B)
    nd <- max(2L, min(n.dist %||% (nb - 1), nb))
    rows <- c(seq_len(nd - 1), nb)
    panel(B[rows, 1], B[rows, -1, drop = FALSE], length(rows), "minimum geodesic distance",
          "proportion of dyads", c(seq_len(nd - 1), "NR"))
  }
  esp_panel <- function() {
    ne <- n.esp %||% length(x$obs.esp)
    panel(x$obs.esp, x$sim.esp, ne, "edge-wise shared partners", "proportion of edges",
          seq(0, ne - 1))
  }
  if (!isTRUE(x$directed)) {
    graphics::par(mfrow = c(1, 3), oma = c(0, 0, 3, 0), mar = c(4, 3, 1.5, 1))
    nd <- n.deg %||% length(x$obs.degree)
    panel(x$obs.degree, x$sim.degree, nd, "degree", "proportion of nodes", seq(0, nd - 1))
    graphics::title("Bayesian goodness-of-fit diagnostics", outer = TRUE)
    dist_panel()
    esp_panel()
  } else {
    graphics::par(mfrow = c(2, 2), oma = c(0, 0, 3, 0), mar = c(4, 3, 1.5, 1))
    ni <- n.ideg %||% length(x$obs.idegree)
    panel(x$obs.idegree, x$sim.idegree, ni, "in degree", "proportion of nodes", seq(0, ni - 1))
    graphics::title("Bayesian goodness-of-fit diagnostics", outer = TRUE)
    no <- n.odeg %||% length(x$obs.odegree)
    panel(x$obs.odegree, x$sim.odegree, no, "out degree", "proportion of nodes", seq(0, no - 1))
    dist_panel()
    esp_panel()
  }
  invisible(x)
}

#' @export
print.bgof <- function(x, ...) {
  cat("Bayesian goodness-of-fit for model: y ~", .formula_rhs(x$formula), "\n")
  cat(ncol(x$sim.esp), "posterior predictive networks.\n\n")
  cat("Model statistics (observed, posterior predictive mean, p-value):\n")
  tab <- cbind(obs = x$obs.model, mean = colMeans(x$sim.model), "p-value" = x$pval.model)
  print(round(tab, 4))
  invisible(x)
}

#' Posterior predictive simulation from a \code{bergm} fit
#'
#' Draws parameters from the posterior sample and simulates one network for
#' each, starting from the observed network.
#'
#' @param object an object of class \code{bergm}.
#' @param nsim count; number of networks.
#' @param seed integer; optional seed.
#' @param aux.iters count; number of MCMC iterations per network.
#' @param output \code{"network"} (a list of networks) or \code{"stats"}
#' (matrix of model statistics).
#' @param cores count; parallel workers (for \code{output = "stats"}).
#' @param ... unused.
#' @export
#' @method simulate bergm
simulate.bergm <- function(object, nsim = 1, seed = NULL, aux.iters = 10000,
                           output = c("network", "stats"), cores = 1, ...) {
  output <- match.arg(output)
  .with_seed(seed)
  Theta <- as.matrix(object$Theta)
  FF <- Theta[sample.int(nrow(Theta), nsim, replace = nsim > nrow(Theta)), , drop = FALSE]
  pp <- .posterior_predictive(object, FF, aux.iters, gof = FALSE, cores = cores,
                              seed = seed, networks = output == "network")
  if (output == "stats") {
    S <- pp$stats
    attr(S, "theta") <- FF
    return(S)
  }
  nets <- pp$networks
  attr(nets, "theta") <- FF
  class(nets) <- "network.list"
  nets
}
