#' Summary of BERGM posterior output
#'
#' This function summarises MCMC output for \code{bergm} objects.
#'
#' @docType methods
#' 
#' @param object an \code{R} object of class \code{bergm}.
#' 
#' @param ... additional arguments, to be passed to lower-level functions.
#' 
#' @export
#' @method summary bergm
 
summary.bergm <- function(object, ...) {
  
  x <- object
  
  stopifnot(inherits(x, "bergm"))

  cat("\n", "Posterior Density Estimate for Model: y ~", paste(x$formula[3]), "\n", "\n")
  
  Theta <- as.mcmc(x$Theta)
  quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)
  
  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
  varstats <- matrix(nrow = nvar(Theta), ncol = length(statnames), 
                     dimnames = list(varnames(Theta), statnames))
  
  Thetavar <- apply(Theta, 2, var)
  Thetatsvar <- apply(Theta, 2, function(x) coda::spectrum0.ar(x)$spec)
  varquant <- t(apply(Theta, 2, quantile, quantiles))
  
  varstats[, 1] <- apply(Theta, 2, mean)
  varstats[, 2] <- sqrt(Thetavar)
  varstats[, 3] <- sqrt(Thetavar / niter(Theta))
  varstats[, 4] <- sqrt(Thetatsvar / niter(Theta))
  
  table1 <- drop(varstats)
  table2 <- drop(varquant)
  
  rNames <- paste("theta", seq(1, x$dim), " (", x$specs[seq(1, x$dim)], ")", sep = "")
  
  rownames(table1) <- rownames(table2) <- rNames
  print(table1); cat("\n"); print(table2)
  cat("\n", "Acceptance rate:", x$AR, "\n", "\n", "\n")
}
