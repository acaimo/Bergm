#' Summary of BERGM posterior output
#'
#' This function summarises MCMC output for \code{bergm} objects.
#'
#' @docType methods
#' 
#' @param object an \code{R} object of class \code{bergm}.
#'
#' @param returnTable logical;
#' By default FALSE. Should the results table be returned.
#' 
#' @param ... additional arguments, to be passed to lower-level functions.
#' 
#' @export
#' @method summary bergm
 

summary.bergm <- function(object,
                          returnTable = FALSE,
                          ...) {
  x <- object

  stopifnot(inherits(x, "bergm"))

  cat("\nPosterior Density Estimate for Model: y ~",
      paste(x$formula[3]),
      "\n\n")

  Theta <- as.mcmc(x$Theta)
  quantiles <- c(0.025, 0.25, 0.5, 0.75, 0.975)

  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE","ESS",
                 "Pr(Theta<0)")
  varstats <- matrix(nrow = nvar(Theta),
                     ncol = length(statnames),
                     dimnames = list(varnames(Theta),
                                     statnames))

  Thetavar <- apply(Theta[,!x$fixed], 2, var)
  Thetatsvar <- coda::spectrum0.ar(Theta[,!x$fixed])$spec
  varquant <- t(apply(Theta[,!x$fixed], 2, quantile, quantiles))

  varstats[, "Mean"] <- apply(Theta, 2, mean)
  varstats[!x$fixed, "SD"] <- sqrt(Thetavar)
  varstats[!x$fixed, "Naive SE"] <- sqrt(Thetavar / niter(Theta))
  varstats[!x$fixed, "Time-series SE"] <- sqrt(Thetatsvar / niter(Theta))
  varstats[, "ESS"] <- x$ess
  varstats[!x$fixed, "Pr(Theta<0)"] <- apply(Theta[,!x$fixed], 2,
                              function(y) {sum(y < 0)/length(y)})

  table1 <- drop(varstats)
  table2 <- drop(varquant)

  print(table1)
  cat("\n")
  print(table2)
  cat("\n", "Acceptance rate:", x$AR, "\n", "\n", "\n")
  if (any(x$fixed)) {
    cat("Some parameters were fixed during the estimation:\n",
        format(row.names(table1)[x$fixed]),'\n\n\n')
  }
  if (returnTable) {
    return(table1)
  }
}
