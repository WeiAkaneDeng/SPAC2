#' A bias corrected criterion for selecting number of principal components
#'
#' The function returns the choice for PCA as a by-product of the bias-corrected residual variance estimate.
#'
#' @param lambda a numeric vector of sample eigenvalues of length $n$.
#' @param M the number of observations.
#' @param constant a small prefixed constant and set to the recommended value of 0.05. See Passemier et al., (2017) for details.
#' @param tau  a tolerance threshold for the smallest eigenvalue, the default value is 0.001.
#'
#' @return an integer K
#'
#' @examples
#' \dontrun{
#' X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' passemier(lambda = eigen_values, M = 100)
#' passemier(lambda = eigen_values, M = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Passemier, D., Li, Z., & Yao, J. (2017). On estimation of the noise variance in high dimensional probabilistic principal component analysis. **Journal of the Royal Statistical Society: Series B (Statistical Methodology)**, *79*(1), 51-67. <doi:10.1111/rssb.12153>
#'
#' @export passemier
#'

passemier <- function(lambda=NULL, M = NULL, tau  = 0.001, constant = 0.05) {

  if (is.null(lambda)) {
    stop("Please provide a numerical vector of sample eigenvalues")
  }


PC_star <- function(k, p, lambda, delta = 0.05, condition=1000, n){
	nn <- length(lambda)
	kmax <- round(nn/2)
	sigma2hat <- function(xx) sum(lambda[(xx+1):nn])/(n-xx)
	c_n <- min(n/p, 1)
	b_sigma2 <- function(xx) sqrt(c_n/2)*(xx + sigma2hat(xx)*sum(1/lambda[1:xx]))
	sigma2hat_star <- function(xx) sigma2hat(xx) + b_sigma2(xx)/(n-xx)*sigma2hat(xx)*sqrt(c_n*2)

sigma2hat_star(k) + sigma2hat_star(kmax)*k*(c_n+2*sqrt(c_n))*(1+p/n^(1+delta))/n
}


	n <- sum(lambda > tau)
	which.min(sapply(1:(n-1), function(x) PC_star(k=x, lambda=lambda, p = M, n = n)))
}



