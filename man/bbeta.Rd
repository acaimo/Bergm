\name{bbeta}
\alias{bbeta}
\title{
Bayesian beta model
}
\description{
Function to fit the Bayesian beta model using an MCMC algorithm.
}
\usage{
bbeta(formula, 
      burn.in = 100, 
      main.iters = 1000, 
      m.prior = NULL, 
      sigma.prior = NULL, 
      sigma.epsilon = NULL,
      ...)	
}
\arguments{
  \item{formula}{
formula; an \code{R} formula object, of the form  <network> ~ <model terms> where <network> is a \code{\link{network}} object and <model terms> are \code{\link{ergm-terms}}.
}
\item{burn.in}{
count; number of burn-in iterations at the beginning of an MCMC run.
}
  \item{main.iters}{
count; number of iterations for the MCMC chain excluding burn-in.
}
  \item{m.prior}{
vector; mean of the multivariate Normal prior. By default set to a vector of 0's.
}
  \item{sigma.prior}{
variance/covariance matrix for the multivariate Normal prior. By default set to a diagonal matrix with every diagonal entry equal to 100.
}
  \item{sigma.epsilon}{
variance/covariance matrix for the multivariate Normal proposal. By default set to a diagonal matrix with every diagonal entry equal to 0.0025.
}
  \item{...}{
additional arguments, to be passed to lower-level functions.
}
}

\examples{
# Make sure the observed network does not have isolated nodes
set.seed(27)

y <- network(7, directed = FALSE)

post.est <- bbeta(y ~ sociality(base = 0), 
                  main.iters = 3000,
                  sigma.epsilon = diag(1.3, dim(y[,])[1]))

bbeta.output(post.est, plot = TRUE)
}