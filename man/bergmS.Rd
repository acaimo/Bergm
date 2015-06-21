\name{bergmS}
\alias{bergmS}
\title{
Bayesian model selection for exponential random graph models
}
\description{
Bayesian model selection for exponential random graphs models using the auto-RJ exchange algorithm. The algorithm consists of two steps: the first step (offline) is used to sample from the posterior of each competing model using the \code{\link{bergm}} function and then approximated by normal distributions determined by the moments of each sample. The second step (online step) of the algorithm makes use of the normal posterior proposal estimated in the offline step as within-model proposals for the RJ-MCMC computation.
}
\usage{
bergmS(formulae, 
       iters = 10000, 
       m.priors = NULL, 
       sigma.priors = NULL, 
       gammas = NULL, 
       nchains = NULL, 
       sigma.epsilons = NULL, 
       aux.iters = 1000, 
       main.iters = NULL,
       burn.ins=NULL,
       ...)
}
\arguments{
  \item{formulae}{
list; \code{R} formula objects (competing models) of the form  <network> ~ <model terms>, see \code{\link{ergm-terms}}.
}
  \item{iters}{
count; number of iterations for auto-RJ exchange algorithm (online step).
}
  \item{m.priors}{
list; vectors of means of the multivariate Normal priors for each competing model. By default set to a list of vectors (one for each competing model) of 0's.
}
  \item{sigma.priors}{
list; variance /covariance matrices of the multivariate Normal priors for each competing model. By default set to a list of diagonal matrices (one for each competing model) with every diagonal entry equal to 100.
}
  \item{gammas}{
vector; ``parallel ADS move factors'' for each competing model (offline step). By default set to a vector of 0.5's.
}
  \item{nchains}{
vector; number of MCMC chains for each competing model (offline step). By default set to a vector with entries (one for each competing model) equal to twice the dimensions of the competing models.
}
  \item{sigma.epsilons}{
list; containing variance/covariance matrices for the multivariate Normal proposals or ``parallel ADS move parameters'' for each competing model (offline step). By default set to a list of diagonal matrices (one for each competing model) with every diagonal entry equal to 0.0025.
}
  \item{aux.iters}{
count; number of auxiliary iterations for network simulation.
}
  \item{main.iters}{
vector; number of iterations for the MCMC chain(s) for each competing model (offline step). By default set to  a vector of 1000's.
}
  \item{burn.ins}{
vector; number of burn-in iterations at the beginning of an MCMC run for each competing model (offline step). By default set to  a vector of 100's.
}
  \item{\dots}{
additional arguments, to be passed to lower-level functions.
}
}
\seealso{
 \code{\link{bergm}}.
}

\examples{
data(florentine)

y <- flomarriage

# Competing models:

formulae <- c(y ~ edges,
              y ~ edges + kstar(2))

# Model selection via auto-RJ approximate exchange algorithm

flo <- bergmS(formulae,
              iters=1500,
              aux.iters=500,
              main.iters=rep(100,2),
              burn.in=rep(50,2),
              gammas=c(1,1))

# MCMC diagnostics

flo.out <- bergmS.output(flo)
}
