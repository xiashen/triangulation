#'
#' Maximum likelihood estimation of sensitivity, specificity, and prevalence
#' 
#' The function estimates sensitivity, specificity, and prevalence parameters, for at least three different methods, without a gold standard. 
#' 
#' @param counts A summary counts vector for the contingencies of each combination of test results, must be named according to the 0-1 pattern, 
#' e.g., 010 = results from 3 tests are negative, positive, and negative, respectively. \code{names(ny)} maybe incomplete and has arbitrary order.
#' @param ntest An integer of the number of tests performed by each method, i.e., sample size for evaluating the likelihood. 
#' @param method The optimization method to be used, can be one of "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", or "Brent". See Note. 
#' e.g., 010 = results from 3 tests are negative, positive, and negative, respectively. \code{names(ny)} maybe incomplete and has arbitrary order.
#' @param B An integer for the number of bootstrap samples, default 0 means no bootstrap estimation for standard errors. 
#' @param start A vector of 2K + 1 parameter starting values, corresponding to K sensitivity values, K specificity values, and prevalence of true positives, respectively. Default: 0.5 for all parameters.
#'
#' @note Please refer to the documentation of the \code{optim} function for method details.
#'
#' 
#' @return A list of parameter estimates. If B > 0, the list also contains the standard errors obtained via bootstrap.
#'
#' 
#' @author Xia Shen, Yudi Pawitan
#' 
#' @references 
#' Yang Z, Xu W, Zhai R, Li T, Ning Z, Pawitan Y, Shen X (2020). Triangulation of analysis strategies links complex traits to specific tissues and cell types. 
#' \emph{Submitted}.
#' 
#' @seealso 
#' \code{combined.fdr}
#'
#' 
#' @examples 
#'\dontrun{
#'
#' ny <- c(1051, 179, 1028, 154, 1040, 159, 981, 208) 
#' names(ny) <- c('000', '001', '010', '011', '100', '101', '110', '111') 
#' ny
#' #  000  001  010  011  100  101  110  111
#' # 1051  179 1028  154 1040  159  981  208
#'
#' triangulate(counts = ny, ntest = 3) # no standard error
#' 
#' triangulate(counts = ny, ntest = 3, B = 10) 
#'
#' }
#' @export
#' 

triangulate <- function(counts, ntest, method = 'Nelder-Mead', B = 0, start = c(rep(.5, ntest*2), .5))
{
	param <- start
	res <- optim(param, loglik, ntest = ntest, ny = counts)$par 
	cat('Estimation done.\n')
	if (B > 0) {
		## bootstrap
		cat('Boostrap standard errors:\n')
		res.boot <- matrix(NA, B, length(param))
		require(svMisc)
		set.seed(911)
		for (j in 1:B) {
			xx <- sample(rep(1:length(counts), counts), sum(counts), replace = TRUE)
			xxtab <- table(xx)
			names(xxtab) <- names(counts)
			res.boot[j,] <- optim(param, loglik, ntest = ntest, ny = xxtab)$par 
			progress(j/B*100)
		}
		se <- sqrt(apply(res.boot, 2, 'var'))
		cat('\n')
		return(list(sensitivity.est = res[1:ntest], sensitivity.se = se[1:ntest], 
		   specificity.est = res[(ntest + 1):(2*ntest)], specificity.se = se[(ntest + 1):(2*ntest)], 
		   prevalence.est = res[2*ntest + 1], prevalence.se = se[2*ntest + 1]))
	} else {
		cat('\n')
		return(list(sensitivity.est = res[1:ntest],  
		   specificity.est = res[(ntest + 1):(2*ntest)], 
		   prevalence.est = res[2*ntest + 1]))
	}	
	
}



