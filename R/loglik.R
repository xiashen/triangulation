#'
#' Log-likelihood of sensitivity, specificity, and prevalence
#' 
#' The function returns the minus log-likelihood value given sensitivity and specificity of at least three testing methods and prevalence of true positives as parameters, 
#' and contingency table of testing results of the methods as data.
#' 
#' @param param A vector of 2K + 1 parameter values, corresponding to K sensitivity values, K specificity values, and prevalence of true positives, respectively.
#' @param ntest An integer of the number of tests performed by each method, i.e., sample size for evaluating the likelihood. 
#' @param ny A summary counts vector for the contingencies of each combination of test results, must be named according to the 0-1 pattern, 
#' e.g., 010 = results from 3 tests are negative, positive, and negative, respectively. \code{names(ny)} maybe incomplete and has arbitrary order.
#' @param tol A value of tolerance, default 1e-4. If any \code{param} value or marginal probability is too close to 0 or 1, distance less than \code{tol}, 
#' a log-likelihood value of \code{exp(30)} will be returned.
#'
#' @note This is an internal function for likelihood calculation.
#' 
#' @return A single value of minus log-likelihood.
#'
#' 
#' @author Xia Shen, Yudi Pawitan
#' 
#' @references 
#' Yang Z, Xu W, Zhai R, Li T, Ning Z, Pawitan Y, Shen X (2020). Triangulation of analysis strategies links complex traits to specific tissues and cell types. 
#' \emph{Submitted}.
#' 
#' @seealso 
#' \code{triangulate}
#'
#' 
#' @examples 
#'\dontrun{
#' sens <- c(.5,.6,.7)
#' spec <- c(.5,.4,.3)
#' p <- .5
#' param <- c(sens, spec, p)
#' ny <- c(123, 456, 789)
#' names(ny) <- c('010', '110', '111')
#' ny
#' # 010 110 111
#' # 123 456 789
#'
#' loglik(param, ntest = 3, ny)
#' # [1] 2844.676
#' }
#' 

loglik <- function(param, ntest, ny, tol = 1e-4) 
{
	if (min(param) <= tol | max(param) >= 1 - tol) {
		loglik <- exp(30);
		return(loglik)
	}
	sens <- param[1:ntest]
	spec <- param[(ntest + 1):(2*ntest)]
	p <- param[2*ntest + 1]
	prob <- mixprob(sens, spec, p) ## this has complete list of names
	##   
	obs.comb <- intersect(names(ny), names(prob))
	obs.ny <- ny[obs.comb]
	model.prob <- prob[obs.comb]
	if (min(prob) > tol & max(prob) < 1 - tol) loglik <- -sum(obs.ny*log(model.prob))
	if (min(prob) <= tol | max(prob) >= 1 - tol) loglik <- exp(30)
	return(loglik)
}
