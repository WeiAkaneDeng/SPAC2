#' A bias corrected criterion for selecting number of principal components
#'
#' The function returns the choice for PCA as a by-product of the bias-corrected residual variance estimate.
#'
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param lambda a numeric vector of sample eigenvalues of the covariance matrix of t(\code{x})
#' @param M if \code{x} were not supplied, \code{M} should be given as the number of columns of \code{x}.
#' @param constant a small prefixed constant and set to the recommended value of 0.05. See Passemier et al., (2017) for details.
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
#' @references Passemier, Damien, Zhaoyuan Li, and Jianfeng Yao (2017). On estimation of the noise variance
#'		in high dimensional probabilistic principal component analysis. \emph{Journal of the Royal Statistical Society:
#'		Series B (Statistical Methodology)} \strong{79.1}: 51-67.
#'
#' @export
#'
passemier <- function(x = NULL, lambda=NULL, M = NULL, constant = 0.05) {

  if (is.null(x) && is.null(lambda)) {
    stop("Please provide either a data matrix or a numerical vector of sample eigenvalues")
  }

  if (is.null(x)) {

    if (is.null(M)) {
      stop("Please provide the number of observations or features along with the sample eigenvalues")
    }

    lambda <- ifelse(lambda > 0, lambda, 0)
    n <- sum(lambda > 0)

    sigma2 <- sapply(1:(n - 1), function(x) sum(lambda[(x + 1):n])/(n - x))

  } else if (is.null(lambda)) {

    X <- x[, !apply(x, 2, function(xx) sum(is.na(xx)) > 0)]
    M <- ncol(X)
    n <- nrow(X)
    S <- stats::cov(scale(t(X)), use = "pairwise.complete")
    SS <- as.matrix(Matrix::nearPD(S)$mat)
    lambda = eigen(SS)$val
    lambda = ifelse(lambda > 0, lambda, 0)
  }

    p <- sum(lambda > 0)
    out <- NA
    for (K in 1:(p-1)){
		sigma2 <- sum(lambda[(K+1):p])/(p-K)
		alpha <- lambda[1:K]
		cn <- p/(M-1)

		b_sigma2 <- sqrt(cn/2)*(K + sigma2*sum(1/alpha[1:K]))
		sigma2_star <- sigma2 + b_sigma2/(p-K)*sigma2*sqrt(2*cn)
		sigma2_star_m <- sigma2 + b_sigma2*sigma2*sqrt(2*cn)
		# max is p-1

		gNT <- (sqrt(cn)+2*sqrt(cn))*(1+M/p^(1+constant))/p
		PC_star <- sigma2_star + K*sigma2_star_m*gNT
		out[K] <- PC_star
		}

      return(min(which.min(out), p-1))
}
