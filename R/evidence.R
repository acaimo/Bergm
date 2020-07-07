#' Wrapper function for evidence estimation
#'
#' Function to estimate the evidence (marginal likelihood) with Chib and Jeliazkov's method
#' or Power posteriors, based on the adjusted pseudolikelihood function.
#' 
#' @param evidence.method vector Method to estimate the marginal likelihood. Options are: \code{"CJ"}, in
#' which case the marginal likelihood is estimated with Chib and Jeliazkov's method; \code{"PP"}, in
#' which case the marginal likelihood is estimated with Power posteriors.
#' 
#' @param ... further arguments to be passed.
#' See \code{evidenceCJ} and \code{evidencePP}.
#'
#' @references
#' Bouranis, L., Friel, N., & Maire, F. (2018). Bayesian model selection for exponential 
#' random graph models via adjusted pseudolikelihoods. 
#' Journal of Computational and Graphical Statistics, 27(3), 516-528. 
#' \url{https://arxiv.org/abs/1706.06344}
#'
#' @examples
#' \dontrun{
#' # Load the florentine marriage network:
#' data(florentine)
#'                                                 
#' # MCMC sampling and evidence estimation:
#' CJE <- evidence(evidence.method = "CJ",
#'                 formula     = flomarriage ~ edges + kstar(2),
#'                 main.iters  = 30000,
#'                 burn.in     = 2000,
#'                 aux.iters   = 1000,
#'                 num.samples = 25000,
#'                 V.proposal  = 2.5,
#'                 ladder      = 100,
#'                 seed        = 1)
#'                                    
#' # Posterior summaries:
#' summary(CJE)
#' 
#' # MCMC diagnostics plots:
#' plot(CJE)
#'     
#' # Log-evidence (marginal likelihood) estimate:
#' CJE$log.evidence
#'}
#'
#' @export
#'
evidence <- function(evidence.method = c("CJ", "PP"),
                        ...)
{
  
  if ( is.element(evidence.method, c("CJ", "PP")) ){
    evidence.method <- match.arg(evidence.method, c("CJ", "PP"))
  } else {
    stop("Select a valid evidence estimation method.\n")
  }
  
  call <- as.list(match.call())[-1]
  
  if( evidence.method == "CJ" ){
    do.call(evidenceCJ, call)
  } else if ( evidence.method == "PP" ){
    do.call(evidencePP, call)
  }
}
