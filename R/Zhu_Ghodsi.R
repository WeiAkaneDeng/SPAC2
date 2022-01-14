#' Automatic Dimensionality Selection from the Scree Plot via Profile Likelihood
#'
#' The function returns the estimated dimension from a scree plot using a profile
#' likelihood criterion. The method requires only the sample eigenvalues.
#'
#' @param lambda a numerical vector of sample eigenvalues of length $n$
#' @param tau  a tolerance threshold for the smallest eigenvalue, the default value is 0.001.
#' @param verbose a logical to indicate whether the details of the profile log-likelihood evidence should be shown
#'
#' @return an integer $K$ between 1 and $n$.
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom mvtnorm dmvnorm
#' @importFrom pracma orth
#' @importFrom psych tr
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' normdata <- mvrnorm(1000, mu = rep(0,50), Sigma = diag(1,50))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(normdata)))$mat))$val
#'
#' ZG(lambda = eigen_values) # supply the sample eigenvalues
#' }
#'
#' @author Wei Q. Deng, \email{dengwq@mcmaster.ca} and Radu V. Craiu \email{craiu@utoronto.ca}
#'
#' @references Zhu, M., & Ghodsi, A. (2006). Automatic dimensionality selection from the scree plot via the use of profile likelihood. **Computational Statistics & Data Analysis**, *51*(2), 918-930. <doi:10.1016/j.csda.2005.09.010>
#'
#' @export ZG


ZG <- function(lambda, tau = 0.001, verbose=FALSE){

  if (is.null(lambda)) {
    stop("Please provide a numerical vector of sample eigenvalues")
  }

  n <- length(lambda > tau)
  prlk <- NA

  for (k in 1:n){

	mu1 <- sum(lambda[1:k])/k;
	mu2 <- sum(lambda[(k+1):n])/(n-k);

	s1 <- stats::var(lambda[1:k]);
	s2 <- stats::var(lambda[(k+1):n]);

	sigma2Hat <- ((k-1)*s1 + (n-k-1)*s2)/(n-2);

	prlk[k] <- sum(stats::dnorm(lambda[1:k], mean = mu1, sd = sqrt(sigma2Hat), log=TRUE)) + sum(stats::dnorm(lambda[(k+1):n], mean = mu2, sd = sqrt(sigma2Hat), log=TRUE))
  }

  if (verbose) {
    return(prlk)
  } else {
    return(which.max(prlk))
  }
	}

