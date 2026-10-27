#' Bayesian exponential random graph models
#'
#' \code{Bergm} provides a range of tools to analyse
#' Bayesian exponential random graph models using advanced
#' computational methods: parameter estimation via the approximate exchange
#' algorithm (\code{\link{bergm}}), calibrated pseudo-posteriors
#' (\code{\link{bergmC}}), adjusted pseudolikelihoods (\code{\link{ergmAPL}}),
#' model evidence estimation (\code{\link{evidence}}), missing data imputation
#' (\code{\link{bergmM}}) and Bayesian goodness-of-fit (\code{\link{bgof}}).
#'
#' @name Bergm-package
#' @aliases Bergm
#' @import ergm
#' @import network
#' @importFrom coda mcmc mcmc.list effectiveSize rejectionRate
#' @importFrom stats coef vcov simulate cov var quantile runif rnorm update
#' @importFrom Rcpp evalCpp
#' @useDynLib Bergm, .registration = TRUE
#' @keywords internal
"_PACKAGE"
