#' Summarising posterior BERGM output
#'
#' This function returns the posterior parameter density estimate 
#' and creates simple diagnostic plots for the MCMC produced from a fit.
#'
#' @param x an \code{R} object of class \code{bergm}, \code{pseudo.bergm}, 
#' or \code{calibrate.bergm}.
#'
#' @param ... additional arguments, to be passed to lower-level functions.
#' 
#' @export
#' 

bergm.output <- function(x, ...) {
  
  if (class(x) == "bergm") {
    if (x$nchains > 1) {
      x$AR <- mean(x$AR)
      FF <- mcmc(apply(x$Theta, 2, cbind))
    } else {
      FF <- mcmc(x$Theta)
      rates <- matrix(x$AR, 1, x$dim)
    }
    
  } else FF <- mcmc(x$Theta)
  
  cat("\n", "Posterior Density Estimate for Model: y ~", paste(x$formula[3]), "\n", "\n")

  table1 <- summary(FF)$statistics
  table2 <- summary(FF)$quantiles
  rownames(table1) <- paste("theta", seq(1, x$dim), 
                            " (", x$specs[seq(1, x$dim)], ")", 
                            sep = "")
  rownames(table2) <- paste("theta", seq(1, x$dim), 
                            " (", x$specs[seq(1, x$dim)], ")", 
                            sep = "")
  
  print(table1); cat("\n"); print(table2)
  cat("\n", "Acceptance rate:", x$AR, "\n", "\n", "\n")
  
  dev.new()
  
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
    plot(density(FF[, i]), 
         main = "", 
         axes = FALSE, 
         xlab = bquote(paste(theta[.(i)], " (", .(x$specs[i]), ")")),
         ylab = "", lwd = 2)
    axis(1)
    axis(2)
    
    traceplot(FF[, i], type = "l", xlab = "Iterations", ylab = "")
    autocorr.plot(FF[, i], auto.layout = FALSE, ...)
    
    if (x$dim > 4) seqq <- seq(4, x$dim, 4)
    if (i %in% union(x$dim, seqq)) {
      title(paste("MCMC output for Model: y ~", x$formula[3]), outer = TRUE)
    }
  }
  out <- list(statistics = summary(FF)$statistics,
              quantiles = summary(FF)$quantiles,
              theta = FF)
}
