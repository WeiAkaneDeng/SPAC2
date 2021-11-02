#' Profile log-likelihood of the PPCA model.
#'
#' The function returns the profile log-likelihood
#'    of the PPCA model for each possible choice of
#'    \eqn{K (1, 2, \dots, n-1)} at their respective MLEs.
#'    The maximum choice was set at \eqn{n-1} because when \eqn{K=n},
#'    the profile log-likelihood is equal to that at \eqn{K=n-1}.
#'
#' @param x a data matrix where the number of rows is to be reduced; only complete columns are used
#' @param lambda a numeric vector of sample eigenvalues of the covariance matrix of t(\code{x})
#' @param M if \code{x} were not supplied, \code{M} should be given as the number of columns of \code{x}.
#'
#' @param EM a logic indicator for whether the profile log-likelihood should be computed for
#'    the MLE estimated using EM algorithms, if FALSE, the profile log-likelihood will be
#'    evaluated by substitute the analytical formulation of MLEs.
#'    Note that if \eqn{n} > \eqn{m}, it will not be possible to use this option.
#'
#' @param param a list of MLEs to be supplied if the profile log-likelihood was to
#'    be evaluate at external MLES; the first list should contain a loading
#'    matrix \eqn{W} with dimension \eqn{n} by \eqn{K} (\eqn{\le n});
#'    and the second a numeric value between 0 and 1 for the
#'    residual variance \code{sigma2}.
#'
#' @return profile log-likelihood of length \eqn{n-1}.
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom mvtnorm dmvnorm
#' @importFrom psych tr
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' U <- pracma::randortho(10)
#' W <- U%*%diag(1, 10)%*%t(U)
#' sigma2 <- 0.4
#' ppcaLog(x = t(X), EM=TRUE) # supply a data matrix
#' ppcaLog(x = t(X), param = list(W, sigma2)) # supply a data matrix and MLEs
#' ppcaLog(lambda = eigen_values, M = 1000) # supply the sample eigenvalues
#' }
#' @author Wei Q. Deng, \email{dengwq@mcmaster.ca}
#'
#' @references Tipping, M. E., and Bishop, C. M. (1999). Probabilistic principal component analysis. **Journal of the Royal Statistical Society: Series B (Statistical Methodology)**, *61*(3), 611-622. <doi:10.1111/1467-9868.00196>
#'
#' @keywords probabilistic PCA, Expectation and Maximization, Maximum Likelihood Estimates, profile log-likelihood
#'
#' @export ppcaLog
#'

ppcaLog <- function(x = NULL, lambda = NULL, M = NULL, param = NULL, EM = FALSE) {

  if (is.null(x) && is.null(lambda)) {
    stop("Please provide either a data matrix or a numerical vector of sample eigenvalues")
  }

  if (is.null(x)) {

    if (is.null(M)) {
      stop("Please provide the number of observations or features along with the sample eigenvalues")
    }

    lambda <- lambda[lambda > 0]
    N <- length(lambda)
    sigma2 <- sapply(1:(N - 1), function(x) sum(lambda[(x + 1):N])/(N - x))
    loglk <- -M/2 * (cumsum(log(lambda)[1:(N - 1)]) + (N - (1:(N - 1))) * log(sigma2) +
                       N + N * log(2 * pi))

  } else if (is.null(lambda)) {

    X <- x[, !apply(x, 2, function(xx) sum(is.na(xx)) > 0)]
    M <- ncol(X)
    N <- nrow(X)
    S <- stats::cov(scale(t(X)), use = "pairwise.complete")
    SS <- as.matrix(Matrix::nearPD(S)$mat)

    if (EM == FALSE){
      if (is.null(param)) {
        stop("Please supply a list of MLEs W and sigma2")
      }
      sigma2M <- param[[2]]
      W <- param[[1]]
      loglk <- -M/2 * (N * log(2 * pi) + log(det(W %*% t(W) + diag(sigma2M, N))) +
                         psych::tr(solve(W %*% t(W) + diag(sigma2M, N)) %*% SS))
    } else{

      evals <- eigen(SS)$val
      evecs <- eigen(SS)$vectors
      lambda <- ifelse(evals > 0, evals, 0)

      loglk <- NA

      for (i in 1:(N - 1)) {
        list_ppca <- ppcaMLE(X, nComp = i, tol = 1e-10)
        W <- list_ppca$W
        sigma2 <- max(list_ppca$sigma2, 0)
        loglk[i] <- -M/2 * (N * log(2 * pi) + log(det(W %*% t(W) + diag(sigma2,N))) +
        psych::tr(solve(W %*% t(W) + diag(sigma2, N)) %*% SS))
      }
    }
  }
  return(loglk)
}



loglk <- function(lam, n, tau=1e-5){
  nn <- sum(lam > tau)
  sigma2_val <- sapply(1:(nn - 1), function(x) sum(lam[(x + 1):nn])/(n - x))
  loglk_original <- -0.5*(cumsum(log(lam[1:(nn-1)])) + (n-(1:(nn-1)))*log(sigma2_val[1:(nn-1)]) + n + n*log(2));
  # only from k = 1 to k = n-2 (standardized)
  loglk_constructed <- c(-0.5*(n+n*log(2)), loglk_original, -0.5*(sum(log(lam[lam > tau])) + n + n*log(2)))
  #plot(loglk_constructed)
  # from k = 0 to k = n-1
  return(list(loglk_original, loglk_constructed))
}

