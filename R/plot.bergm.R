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
#' @export
#' @method plot bergm
 
plot.bergm <- function(x, ...) {
  
  stopifnot(inherits(x, "bergm"))
  
  seqq <- 4
  par(mfrow = c(min(4, x$dim), 3), 
      oma   = c(0, 0, 3, 0), 
      mar   = c(4, 3, 1.5, 1))
  
  for (i in 1:x$dim) {
    if (i %in% c(5, 9, 13)) {
      dev.new()
      par(mfrow = c(min(4, x$dim - (i - 1)), 3), 
          oma = c(0, 0, 3, 0), 
          mar = c(4, 3, 1.5, 1))
    }
    plot(density(x$Theta[, i]), 
         main = "", 
         axes = FALSE, 
         xlab = bquote(paste(theta[.(i)], " (", .(x$specs[i]), ")")),
         ylab = "", lwd = 2)
    axis(1)
    axis(2)
    traceplot(x$Theta[, i], type = "l", xlab = "Iterations", ylab = "")
    autocorr.plot(x$Theta[, i], auto.layout = FALSE, ...)
    if (x$dim > 4) seqq <- seq(4, x$dim, 4)
    if (i %in% union(x$dim, seqq)) {
      title(paste("MCMC output for Model: y ~", x$formula[3]), outer = TRUE)
    }
  }
}
