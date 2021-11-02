#' Bai and Ng's Criteria for Selecting Number of Probabilistic Factors
#'
#' The function returns the choice for the number of factors in the context of probabilistic factor analysis.
#'
#' @param lambda a numeric vector of sample eigenvalues of length $n$.
#' @param M the number of observations.
#' @param tau a tolerance threshold for the smallest eigenvalue, the default value is 0.001.
#' @param option an integer specifying the choice of the preferred asymptotics regime. See Bai and Ng (2002) for details.
#' @param bias a logical specifying whether the residual variance estimate used in the criterion should be biased or unbiased. See Bai and Ng (2002) for details.
#' @param verbose a logical specifying whether the posterior evidence or
#'     the integer that minimized the evidence should be returned
#'
#' @return an integer K between 1 and n
#'
#' @examples
#' \dontrun{
#' X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' BaiNg(lambda = eigen_values, M=100)
#' BaiNg(lambda = eigen_values, M=5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Bai, J., & Ng, S. (2002). Determining the number of factors in approximate factor models. **Econometrica**, *70*(1), 191-221. <doi:10.1111/1468-0262.00273>
#'
#' @export BaiNg
#'

BaiNg <- function(lambda, M, bias = TRUE, option = 1, tau = 0.001, verbose = FALSE){


Bai <- function(k, p, lambda, bias = TRUE, option = 1, n){

  nn <- length(lambda > tau)
  kmax = round(nn/2)
	sigma2hat <- function(xx) sum(lambda[(xx+1):nn])/(n-xx)
	c_n <- n/p
	b_sigma2 <- function(xx) sqrt(c_n/2)*(xx + sigma2hat(xx)*sum(1/lambda[1:xx]))
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

  	n <- length(lambda > tau)
  	BaiNgs <- sapply(1:(n-1), function(x) Bai(k=x, lambda=lambda, p=M, bias = bias, option = option, n = length(lambda)));

  	if (verbose) {
  	  return(BaiNgs)
  	} else {
  	  return(which.max(BaiNgs))
  	}
}

