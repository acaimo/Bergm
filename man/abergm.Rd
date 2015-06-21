\name{abergm}
\alias{abergm}
\title{
Adaptive algorithms for Bayesian exponential random graph models
}
\description{
Function to fit Bayesian exponential random graphs models using adaptive exchange algorithms.}
\usage{
abergm(formula, 
      burn.in = 100, 
      main.iters = 1000, 
      aux.iters = 1000, 
      m.prior = NULL, 
      sigma.prior = NULL, 
      nchains = NULL, 
      gamma = 0.5, 
      method = 'Adaptive.chains', 
      rectangular = TRUE,
      sigma.epsilon = NULL, 
      updategap = 10,
      ...)
}
\arguments{
  \item{formula}{
formula; an \code{R} formula object, of the form  <network> ~ <model terms> where <network> is a \code{\link{network}} object and <model terms> are \code{\link{ergm-terms}}.
}
\item{burn.in}{
count; number of burn-in iterations at the beginning of an MCMC run. If population MCMC is performed, it refers to the number of burn-in iterations for every chain of the population.
}
  \item{main.iters}{
count; number of iterations for the MCMC chain(s) excluding burn-in. If population MCMC is performed, it refers to the number of iterations for every chain of the population.
}
  \item{aux.iters}{
count; number of auxiliary iterations used for network simulation.
}
  \item{m.prior}{
vector; mean of the multivariate Normal prior. By default set to a vector of 0's.
}
  \item{sigma.prior}{
variance/covariance matrix for the multivariate Normal prior. By default set to a diagonal matrix with every diagonal entry equal to 100.
}
  \item{nchains}{
count; number of chains of the population MCMC. By default set to twice the model dimension (number of model terms). If the model is one-dimensional, \code{nchains} is set to 1.
}
  \item{gamma}{
scalar; ``parallel ADS move factor.'' In case of one-dimensional models, the population MCMC procedure is disabled and \code{gamma} is used as variance of the Normal proposal distribution.
}
  \item{method}{
name of the adaptive strategy: \code{ADS} = adaptive direction sampling, \code{Adaptive.past} = adaptive strategy where past parameter particles are used, \code{Adaptive.chains} (default) = adaptive strategy all particles at the current time for all chains are used.
}
  \item{rectangular}{
logical; if \code{Adaptive.past} is used, it defines the type of adaptive strategy: if \code{TRUE} (default) = all parameter particles from all chains and all past simulations are used, if \code{FALSE} = all parameter particles along the same chain and all are used.
}
  \item{sigma.epsilon}{
variance/covariance matrix for the multivariate Normal proposal or ``parallel ADS move parameter''. By default set to a diagonal matrix with every diagonal entry equal to 0.0025. If the model is one-dimensional, \code{sigma.espilon} is set equal to \code{gamma}.
}
  \item{updategap}{
scalar; iteration interval used to update the variance/covariance matrix.
}
  \item{...}{
additional arguments, to be passed to lower-level functions.
}
}
\seealso{
 \code{\link{bergm}}
}

\examples{
data(molecule)

mol <- abergm(molecule ~ edges + kstar(2),
              burn.in = 50,
              aux.iters = 50,
              main.iters = 500,
              method = 'Adaptive.chains',
              nchains = 4,
              gamma = 1)

}             
              
