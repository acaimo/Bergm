#' Summary of BERGM posterior output
#'
#' This function summarises MCMC output for \code{bergm} objects.
#'
#' @param object an \code{R} object of class \code{bergm}.
#'
#' @param returnTable logical; by default \code{FALSE}. Should the results
#' table be returned.
#'
#' @param ... additional arguments, to be passed to lower-level functions.
#'
#' @return Invisibly (or, with \code{returnTable = TRUE}, visibly) the table of
#' posterior summaries. When per-chain samples are available, the table also
#' reports the potential scale reduction factor (\code{Rhat}).
#'
#' @export
#' @method summary bergm
summary.bergm <- function(object, returnTable = FALSE, ...) {
  x <- object
  stopifnot(inherits(x, "bergm"))

  cat("\nPosterior Density Estimate for Model: y ~", .formula_rhs(x$formula), "\n\n")

  Theta <- as.matrix(x$Theta)
  fixed <- .fixed_of(x)
  free <- !fixed
  quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE", "ESS", "Pr(Theta<0)")
  has_rhat <- !is.null(x$Rhat) && any(is.finite(x$Rhat))
  if (has_rhat) statnames <- c(statnames, "Rhat")
  varstats <- matrix(NA_real_, ncol(Theta), length(statnames),
                     dimnames = list(colnames(Theta) %||% x$specs, statnames))
  n <- nrow(Theta)
  Tf <- Theta[, free, drop = FALSE]
  v <- apply(Tf, 2, stats::var)
  tsv <- if (is.null(x$chains)) coda::spectrum0.ar(Tf)$spec
  varstats[, "Mean"] <- colMeans(Theta)
  varstats[free, "SD"] <- sqrt(v)
  varstats[free, "Naive SE"] <- sqrt(v / n)
  varstats[free, "Time-series SE"] <- if (!is.null(x$chains) && !is.null(x$ess)) {
    sqrt(v / pmax(x$ess[free], 1))
  } else sqrt(tsv / n)
  varstats[, "ESS"] <- x$ess %||% NA
  varstats[free, "Pr(Theta<0)"] <- colMeans(Tf < 0)
  if (has_rhat) varstats[, "Rhat"] <- round(x$Rhat, 3)
  varquant <- t(apply(Tf, 2, stats::quantile, quantiles))

  print(varstats)
  cat("\n")
  print(varquant)
  cat("\n", "Acceptance rate:", x$AR, "\n", "\n", "\n")
  if (any(fixed)) {
    cat("Some parameters were fixed during the estimation:\n",
        format(rownames(varstats)[fixed]), "\n\n\n")
  }
  if (!is.null(x$log.evidence)) cat(" Log-evidence estimate:", format(x$log.evidence), "\n\n")
  if (returnTable) return(varstats)
  invisible(varstats)
}

#' @export
print.bergm <- function(x, ...) {
  cat("Bayesian ERGM fit: y ~", .formula_rhs(x$formula), "\n")
  fixed <- .fixed_of(x)
  m <- colMeans(as.matrix(x$Theta))
  sds <- apply(as.matrix(x$Theta), 2, stats::sd)
  sds[fixed] <- NA
  tab <- cbind(Mean = m, SD = sds)
  rownames(tab) <- x$specs
  print(round(tab, 4))
  cat(sprintf("\n%d draws%s; acceptance rate %s; runtime %s\n", nrow(x$Theta),
              if (!is.null(x$nchains)) sprintf(" from %d chains", x$nchains) else "",
              format(x$AR), format(x$Time, digits = 3)))
  if (!is.null(x$log.evidence)) cat("Log-evidence:", format(x$log.evidence), "\n")
  invisible(x)
}

.fixed_of <- function(x) {
  f <- x$fixed
  if (is.null(f)) f <- rep(FALSE, ncol(as.matrix(x$Theta)))
  as.logical(f)
}

#' Plot BERGM posterior output
#'
#' This function creates MCMC diagnostic plots for \code{bergm} objects:
#' posterior density, trace plot (one line per chain when available) and
#' autocorrelation for every free parameter.
#'
#' @param x an \code{R} object of class \code{bergm}.
#'
#' @param lag.max integer; maximum lag for the autocorrelation plot.
#'
#' @param ... additional arguments, to be passed to lower-level functions.
#'
#' @examples
#' \dontrun{
#' data(florentine)
#' p.flo <- bergm(flomarriage ~ edges + kstar(2), main.iters = 1000, seed = 1)
#' plot(p.flo)
#' }
#'
#' @export
#' @method plot bergm
plot.bergm <- function(x, lag.max = 50, ...) {
  stopifnot(inherits(x, "bergm"))
  free <- which(!.fixed_of(x))
  Theta <- as.matrix(x$Theta)
  specs <- x$specs %||% colnames(Theta)
  real_dim <- length(free)
  per_page <- 4L
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  cols <- if (!is.null(x$chains)) grDevices::gray.colors(length(x$chains), start = 0.1, end = 0.7) else "black"
  for (start in seq(1L, real_dim, by = per_page)) {
    idx <- free[start:min(real_dim, start + per_page - 1L)]
    if (start > 1L && grDevices::dev.interactive()) grDevices::dev.new()
    graphics::par(mfrow = c(length(idx), 3), oma = c(0, 0, 3, 0), mar = c(4, 3, 1.5, 1))
    for (i in idx) {
      th <- Theta[, i]
      graphics::plot(stats::density(th), main = "", axes = FALSE,
                     xlab = bquote(paste(theta[.(i)], " (", .(specs[i]), ")")),
                     ylab = "", lwd = 2)
      graphics::axis(1)
      graphics::axis(2)
      if (!is.null(x$chains)) {
        M <- sapply(x$chains, function(ch) as.matrix(ch)[, i])
        graphics::matplot(M, type = "l", lty = 1, col = cols, xlab = "Iterations", ylab = "")
      } else {
        graphics::plot(th, type = "l", xlab = "Iterations", ylab = "")
      }
      series <- if (!is.null(x$chains)) lapply(x$chains, function(ch) as.matrix(ch)[, i]) else list(th)
      acs <- sapply(series, function(v) {
        lm <- min(lag.max, length(v) - 1L)
        a <- if (stats::var(v) > 0) drop(stats::acf(v, lag.max = lm, plot = FALSE)$acf) else c(1, rep(0, lm))
        c(a, rep(NA, lag.max - lm))
      })
      acs <- rowMeans(matrix(acs, nrow = lag.max + 1L), na.rm = TRUE)
      graphics::plot(seq_along(acs) - 1L, acs, type = "h", ylim = c(-1, 1), xlab = "Lag",
                     ylab = "Autocorrelation", ...)
      graphics::abline(h = 0)
    }
    graphics::title(paste("MCMC output for Model: y ~", .formula_rhs(x$formula)), outer = TRUE)
  }
  invisible(x)
}

#' Posterior point estimates and covariance for \code{bergm} objects
#'
#' @param object an object of class \code{bergm}.
#' @param ... unused.
#' @return \code{coef}: posterior means; \code{vcov}: posterior covariance matrix.
#' @name bergm-methods
NULL

#' @rdname bergm-methods
#' @export
coef.bergm <- function(object, ...) {
  m <- colMeans(as.matrix(object$Theta))
  names(m) <- object$specs
  m
}

#' @rdname bergm-methods
#' @export
vcov.bergm <- function(object, ...) {
  V <- stats::cov(as.matrix(object$Theta))
  dimnames(V) <- list(object$specs, object$specs)
  V
}

#' @rdname bergm-methods
#' @param x an object of class \code{bergm}.
#' @importFrom coda as.mcmc
#' @export
as.mcmc.bergm <- function(x, ...) x$Theta

#' @rdname bergm-methods
#' @importFrom coda as.mcmc.list
#' @method as.mcmc.list bergm
#' @export
as.mcmc.list.bergm <- function(x, ...) {
  if (!is.null(x$chains)) x$chains else coda::mcmc.list(x$Theta)
}
