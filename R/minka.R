#' Automatic Choice of dimensionality for PCA using Laplace's method
#'
#' The function returns the choice dimension for PCA under the PPCA setup using Laplace approximation.
#'
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param lambda a numeric vector of sample eigenvalues of the covariance matrix of t(\code{x})
#' @param M if \code{x} were not supplied, \code{M} should be given as the number of columns of \code{x}.
#' @param evidence a logical specifying whether the BIC evidence or
#'     the integer that minimized the evidence should be returned
#' @return an integer K between 1 and N that maximizes the posterior BIC approximated by Laplace's method.
#'
#' @examples
#' \dontrun{
#' X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' minka2001(lambda = eigen_values, p = 100)
#' minka2001(lambda = eigen_values, p = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Minka, T. P. (2000). Automatic choice of dimensionality for PCA.
#' In \emph{NIPS} (Vol. \strong{13}, pp. 598-604).
#'
#' @export
#'
minka2001 <- function(x = NULL, lambda=NULL, p = NULL, evidence = FALSE) {

  if (is.null(lambda)) {
    stop("Please provide either a data matrix or a numerical vector of sample eigenvalues")
  }


   if (is.null(p)) {
      stop("Please provide the number of observations or features along with the sample eigenvalues")
    }

    lambda <- ifelse(lambda > 0, lambda, 0)
    n <- sum(lambda > 0)

    sigma2 <- sapply(1:(n - 1), function(x) sum(lambda[(x + 1):n])/(n - x))



logDk <- function(k, lambda, p){

	lambda <- ifelse(lambda > 1e-5, lambda, NA)
	n <- length(lambda)
	sigma2 <- sum(lambda[(k+1):n], na.rm=T)/(n-k)

	small_sumN <- NA
	for (ii in 1:k){
		small_sumN[ii] <- sum(log(lambda[ii]-lambda[(ii+1):n]), na.rm=T) + (n-ii)*log(1/sigma2 - 1/lambda[ii])
		}

-p/2*sum(log(lambda[1:k])) - p*(n-k)/2*log(sigma2) - (sum(small_sumN) + log(p)*(2*n*k-k^2-k)/2 )/2 -k/2*log(p) + sum(lgamma((n-(1:k)+1)/2)) +  (2*n*k-k^2-3*k)/4*log(2) + 3*k/2*log(2)

}

 
 	n <- length(lambda)
	Minka <- sapply(1:(n-1), function(x) logDk(k=x, lambda= lambda, p = p)))



    if (evidence) {
      return(Minka)
    } else {
      return(which.max(Minka))
    }
}





#' Automatic Choice of dimensionality for PCA using BIC approximation.
#'
#' The function returns the choice dimension for PCA under the PPCA setup using a simplification of Laplace's method.
#'
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param lambda a numeric vector of sample eigenvalues of the covariance matrix of t(\code{x})
#' @param M if \code{x} were not supplied, \code{M} should be given as the number of columns of \code{x}.
#' @param evidence a logical specifying whether the BIC evidence or
#'     the integer that minimized the evidence should be returned
#' @return an integer K between 1 and N that maximizes the approximated marginal log-likelihood.
#'
#' @examples
#' \dontrun{
#' X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale((X))))$mat))$val
#' minka2001_BIC(lambda = eigen_values, p = 1000)
#' minka2001_BIC(lambda = eigen_values, p = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Minka, T. P. (2000). Automatic choice of dimensionality for PCA.
#' In \emph{NIPS} (Vol. \strong{13}, pp. 598-604).
#'
#' @export

minka2001_BIC <- function(lambda=NULL, p = NULL,evidence = FALSE) {

  if (is.null(lambda)) {
    stop("Please provide  a numerical vector of sample eigenvalues")
  }

 BIClog <- function(lambda, p){
 	 	n <- length(lambda)

 	-2*profilelog(lambda= lambda, p = p) + log(p)*((1:(n-1))*n+(1:(n-1))/2-(1:(n-1))^2/2)/2

 }


	bic <- BIClog(lambda,p=p)

     if (evidence) {
      return(bic)
    } else {
      return(which.min(bic))
    }
}



profilelog <- function(lambda, p){

 	n <- length(lambda)


 	sapply(1:(n-1), function(k) -p/2*(sum(log(lambda)[1:k]) + (n-k)*log(sum(lambda[(k + 1):n])/(n - k)) + n*log(2*pi) + n))

 }


