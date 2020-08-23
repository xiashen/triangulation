#'
#' Calculating combined false positive rate based on Monte Carlo sampling
#' 
#' The function returns the estimated false discovery rate for each category of testing results, given the operating characteristics of different testing methods.
#' 
#' @param sens A vector of K sensitivity values for K different testing methods.
#' @param specs A vector of K specificity values for K different testing methods, with the order of the methods same as in \code{sens}.
#' @param prev Prevalence of true positives.
#' @param size A large integer for Monte Carlo sampling. Default: 1 million.
#'
#' @note The function is based on Monte Carlo sampling estimation. Make sure a sufficiently large \code{size} is given.
#' 
#' @return FDR estimates for each category of testing results combination. 
#'
#' 
#' @author Zhijian Yang, Xia Shen
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
#'
#' combined.fdr(c(.5,.6,.7), c(.5,.4,.3), .5)
#' #        000        001        010        011        100        101        110        111
#' # 0.50135036 0.50012141 0.49931784 0.49982606 0.49945863 0.49953203 0.50014472 0.50040443
#' }
#' @export
#' 

combined.fdr <- function(sens, specs, prev, size = 1e6) 
{
	P <- round(prev*size)
	N <- round((1 - prev)*size)
	real <- c(rep(0, N), rep(1, P))

	result <- real
	for (i in 1:length(sens)) {
		sen <- sens[i]
		spec <- specs[i]
		tn <- round(spec*N)
		fp <- round((1 - spec)*N)
		tp <- round(sen*P)
		fn <- round((1 - sen)*P)
		# N = tn + fp
		method.N <- c(rep(0, tn), rep(1, fp))
		# P = tp + fn
		method.P <- c(rep(1, tp), rep(0, fn))
		method.N <- sample(method.N)
		method.P <- sample(method.P)
		method <- c(method.N, method.P)
		result <- cbind(result, method)
		colnames(result)[i + 1] <- paste0("method", i)
		#fdr <- (1 - spec) * (1 - Prev)/((1 - spec) * (1-Prev) + sen * Prev)
		# cat(paste0("method",i,"\n"))
		# cat(paste0("fdr:",fdr))
		# print(table(result[,1],result[,i+1]))
	}
	colnames(result)[1] <- "real" 
	# print(head(result))

	combination <- result[,2]
	for (i in 3:ncol(result)) {
		combination <- paste0(combination, result[,i])
	}

	tab <- table(result[,'real'], combination)
	# print(tab)
	combined_FDR <- c()
	for (i in 1:ncol(tab)) {
		combined_FDR[i] <- tab[1,i]/(tab[1,i] + tab[2,i])
	}
	names(combined_FDR) <- colnames(tab)
	return(combined_FDR)
}

