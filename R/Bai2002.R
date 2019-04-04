#' Bai and Ng's criteria for selecting number of probabilistic factors
#'
#' The function returns the choice for the number of factors in the context of probabilistic factor analysis.
#'
#' @param lambda a numeric vector of sample eigenvalues of length $n$.
#' @param p the number of observations to calculate the sample eigenvalues.
#' @param choice an integer specifying the choice of the preferred asymptotics. See Bai and Ng (2002s) for details.
#' @param unbiased a logical specifying whether the residual variance estimate used in the criterion should be biased or unbiased.
#'    See Bai and Ng (2002) for details.
#'
#' @return an integer K
#'
#' @examples
#' \dontrun{
#' X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' Bai2002(lambda = eigen_values, p = 100)
#' Bai2002(lambda = eigen_values, p = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Bai, Jushan, and Serena Ng (2002). Determining the number of factors in approximate factor models.
#'    \emph{Econometrica} \strong{70.1}: 191-221.
#'
#' @export
#'

BaiNg <- function(lambda, p, bias = T, option = 1){
	
	Fp <- function(k, p, lambda, bias = T, option = 1){
	
	n <- length(lambda)
	kmax <- round(n/2)
	sigma2hat <- function(xx) sum(lambda[(xx+1):n])/(n-xx)
	c_n <- n/p
	b_sigma2 <- function(xx) sqrt(c_n/2)*(xx+ sigma2hat(xx)*sum(1/lambda[1:xx]))
	sigma2hat_star <- function(xx) sigma2hat(xx) + b_sigma2(xx)/(n-xx)*sigma2hat(xx)*sqrt(c_n*2)
	
	sigma_est <- function(xx) ifelse(bias, sigma2hat_star(xx), sigma2hat(xx))
	
	if (option == 1) {
		
		sigma_est(k) + sigma_est(kmax)*k*(p+n)/n/p*log(n*p/(n+p))
	
	} else if (option == 2) {
		
		sigma_est(k) + sigma_est(kmax)*k*(p+n)/n/p*log(min(n, p))

	} else {
		
		sigma_est(k) + sigma_est(kmax)*k*log(min(n, p))/min(n, p)

	}	
}


	n <- length(lambda)
  	BaiNgs <- which.min(sapply(1:(n-1), function(x) Bai(k=x, lambda=lambda, p = p, bias = bias, option = option)));

	return(BaiNgs)
	
}