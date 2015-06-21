\name{bergm}
\alias{bergm}
\title{
Bayesian exponential random graph models
}
\description{
Function to fit Bayesian exponential random graphs models using the exchange algorithm. Two are the sampling approaches available: block update and population MCMC with parallel Adaptive Direction Sampling (ADS).
}
\usage{
bergm(formula, 
      burn.in=10,
      main.iters=1000,
      aux.iters=1000, 
      m.prior = NULL, 
      sigma.prior = NULL, 
      nchains = NULL, 
      gamma = 0.5, 
      sigma.epsilon = NULL,
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
  \item{sigma.epsilon}{
variance/covariance matrix for the multivariate Normal proposal or ``parallel ADS move parameter''. By default set to a diagonal matrix with every diagonal entry equal to 0.0025. If the model is one-dimensional, \code{sigma.espilon} is set equal to \code{gamma}.
}
  \item{...}{
additional arguments, to be passed to lower-level functions.
}
}
\seealso{
 \code{\link{bergm.output}}, 
 \code{\link{bgof}}.
}

\examples{
data(florentine)

# Parameter estimation via approximate exchange algorithm with ADS:

flo <- bergm(flomarriage ~ edges + kstar(2),
             burn.in=10,
             aux.iters=500,
             main.iters=500,
             gamma=1)

# MCMC diagnostics:

bergm.output(flo)

# Bayesian goodness-of-fit test:

bgof(flo,
     aux.iters=500,
     sample.size=50,
     n.deg=10,
     n.dist=9,
     n.esp=6)
}