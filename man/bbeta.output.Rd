\name{bbeta.output}
\alias{bbeta.output}
\title{
MCMC output diagnostics for the Bayesian beta model estimation procedure
}
\description{
This function returns the posterior parameter density estimate and creates simple diagnostic plots for the MCMC produced from a fit. 
}
\usage{
bbeta.output(x, 
             plot = FALSE,
             ...)
}
\arguments{
  \item{x}{ 
  an \code{R} object of class \code{bbeta}.
}
  \item{plot}{ 
  logical; if \code{TRUE} diagnostic plots for each parameter are plotted.
}
  \item{...}{
additional arguments, to be passed to lower-level functions.
}
}
\seealso{
  \code{\link{bbeta}}.
}
