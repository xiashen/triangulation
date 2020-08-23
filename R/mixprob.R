#'
#' Marginal probability as a function of sensitivity, specificity, and prevalence
#' 
#' The function returns the marginal probability values based on sensitivity and specificity of at least three testing methods and prevalence of true positives.
#' 
#' @param sens A vector of K sensitivity values for K different testing methods.
#' @param spec A vector of K specificity values for K different testing methods, with the order of the methods same as in \code{sens}.
#' @param p Prevalence of true positives.
#'
#' @note This is an internal function for likelihood calculation.
#' 
#' @return A named vector of probabilities for the categories of testing results. The naming is from '00...0' to '11...1', with length K for the K methods, 
#' 0-1 pattern corresponds to the test results, e.g., 010 = results from 3 tests are negative, positive, and negative, respectively. 
#' The order is the same as in input \code{sens} and \code{spec}. 
#'
#' 
#' @author Xia Shen, Yudi Pawitan
#' 
#' @references 
#' Yang Z, Xu W, Zhai R, Li T, Ning Z, Pawitan Y, Shen X (2020). Triangulation of analysis strategies links complex traits to specific tissues and cell types. 
#' \emph{Submitted}.
#' 
#' @seealso 
#' \code{loglik}
#'
#' 
#' @examples 
#'\dontrun{
#' mixprob(sens = c(.5,.6,.7), spec = c(.5,.4,.3), p = .5)
#'
#' #  000  100  010  110  001  101  011  111
#' # 0.06 0.06 0.09 0.09 0.14 0.14 0.21 0.21
#' }
#' 

mixprob <- function(sens, spec, p) 
{
	ntest <- length(sens)
	a <- list()
	for (i in 1:ntest) a[[i]] <- 0:1
	ytab <- expand.grid(a)
	cname <- apply(ytab, 1, paste0, collapse = '')
	##
	prob1 <- t(sens**t(ytab)*(1 - sens)**t(1 - ytab))
	prob1 <- apply(prob1, 1, prod)
	##
	prob0 <- t((1 - spec)**t(ytab)*spec**t(1 - ytab))
	prob0 <- apply(prob0, 1, prod)
	mixprob = p*prob1 + (1 - p)*prob0 
	names(mixprob) <- cname
	return(mixprob)
}
