#' Plot BERGM posterior output
#'
#' This function creates  MCMC diagnostic plots for \code{bergm} objects.
#'
#' @docType methods
#'
#' @param x an \code{R} object of class \code{bergm}.
#'
#' @param ... additional arguments, to be passed to lower-level functions.
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
#' # MCMC diagnostics plots:
#' plot(p.flo)
#' }
#'
#' @export
#' @method plot bergm

plot.bergm <- function(x, ...) {

  stopifnot(inherits(x, "bergm"))
  real_dim <- sum(!x$fixed)
  seqq <- 4
  par(mfrow = c(min(4, real_dim), 3),
      oma   = c(0, 0, 3, 0),
      mar   = c(4, 3, 1.5, 1))

  for (i in 1:real_dim) {
    if (i %in% c(5, 9, 13)) {
      dev.new()
      par(mfrow = c(min(4, real_dim - (i - 1)), 3),
          oma = c(0, 0, 3, 0),
          mar = c(4, 3, 1.5, 1))
    }
    plot(density(x$Theta[,!x$fixed][,i]),
         main = "",
         axes = FALSE,
         xlab = bquote(paste(theta[.(i)], " (", .(x$specs[!x$fixed][i]), ")")),
         ylab = "", lwd = 2)
    axis(1)
    axis(2)
    traceplot(x$Theta[,!x$fixed][,i], type = "l", xlab = "Iterations", ylab = "")
    autocorr.plot(x$Theta[,!x$fixed][,i], auto.layout = FALSE, ...)
    if (real_dim > 4) seqq <- seq(4, real_dim, 4)
    if (i %in% union(real_dim, seqq)) {
      title(paste("MCMC output for Model: y ~", x$formula[3]), outer = TRUE)
    }
  }
}
