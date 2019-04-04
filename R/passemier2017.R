#' A bias corrected criterion for selecting number of principal components
#'
#' The function returns the choice for PCA as a by-product of the bias-corrected residual variance estimate.
#'
#' @param lambda a numeric vector of sample eigenvalues of length $n$.
#' @param p the number of observations to calculate the sample eigenvalues.
#' @param constant a small prefixed constant and set to the recommended value of 0.05. See Passemier et al., (2017) for details.
#'
#' @return an integer K
#'
#' @examples
#' \dontrun{
#' X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' passemier(lambda = eigen_values, p = 100)
#' passemier(lambda = eigen_values, p = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Passemier, Damien, Zhaoyuan Li, and Jianfeng Yao (2017). On estimation of the noise variance
#'		in high dimensional probabilistic principal component analysis. \emph{Journal of the Royal Statistical Society:
#'		Series B (Statistical Methodology)} \strong{79.1}: 51-67.
#'
#' @export
#'
passemier <- function(lambda=NULL, p = NULL, constant = 0.05) {

  if (is.null(lambda)) {
    stop("Please provide a numerical vector of sample eigenvalues")
  }

	PC_star <- function(k, p, lambda, delta = 0.05){
	
	n <- length(lambda)
	kmax <- round(n/2)
	sigma2hat <- function(xx) sum(lambda[(xx+1):n])/(n-xx)
	c_n <- n/p
	b_sigma2 <- function(xx) sqrt(c_n/2)*(xx+ sigma2hat(xx)*sum(1/lambda[1:xx]))
	sigma2hat_star <- function(xx) sigma2hat(xx) + b_sigma2(xx)/(n-xx)*sigma2hat(xx)*sqrt(c_n*2)
		
	sigma2hat_star(k) + sigma2hat_star(kmax)*k*(c_n+2*sqrt(c_n))*(1+p/n^(1+delta))/n
}

	n <- length(lambda)
	which.min(sapply(1:(n-1), function(x) PC_star(k=x, lambda=lambda, p = p)))
}



