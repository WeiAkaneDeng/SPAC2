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
#' minka2001(lambda = eigen_values, M = 100)
#' minka2001(lambda = eigen_values, M = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Minka, T. P. (2000). Automatic choice of dimensionality for PCA.
#' In \emph{NIPS} (Vol. \strong{13}, pp. 598-604).
#'
#' @export
#'
minka2001 <- function(x = NULL, lambda=NULL, M = NULL, evidence = FALSE) {

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

    N <- sum(lambda > 0)
    out <- NA
    for (K in 1:N){
    mm <- N * K - (K + 1) * K/2
    sigma2 <- sum(lambda[(K + 1):N])/(N - K)
    l_tild <- ifelse(lambda > sigma2, lambda, sigma2)
    tiny_prod <- function(i) {
        prod(1/l_tild[(i + 1):N] - 1/l_tild[i]) * prod(lambda[i] - lambda[(i + 1):N])
    }
    logA <- mm * log(M) + sum(log(sapply(1:K, tiny_prod)))
    logrho_k <- 1/(N - K) * sum(log(lambda[(K + 1):N])) - 1/(N - K) * log(sum(lambda[(K +
        1):N]))
    out[K] <- (mm - K)/2 * log(2) - M * (N - K) * logrho_k - 0.5 * logA - K/2 * log(M) + sum(lgamma((N -
        (1:K) + 1)/2)) - M/2 * sum(log(lambda[1:K])) + M * sum(log(lambda[(K + 1):N]))
    }
    if (evidence) {
      return(out)
    } else {
      return(which.max(out))
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
#' minka2001_BIC(lambda = eigen_values, M = 1000)
#' minka2001_BIC(lambda = eigen_values, M = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Minka, T. P. (2000). Automatic choice of dimensionality for PCA.
#' In \emph{NIPS} (Vol. \strong{13}, pp. 598-604).
#'
#' @export

minka2001_BIC <- function(x = NULL, lambda=NULL, M = NULL,evidence = FALSE) {

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

    out <- NA
    N <- sum(lambda > 0)
    for (K in 1:N){
    mm <- N * K - (K + 1) * K/2
    sigma2 <- sum(lambda[(K + 1):N])/(N - K)
    l_tild <- ifelse(lambda > sigma2, lambda, sigma2)
    tiny_prod <- function(i) {
        prod(1/l_tild[(i + 1):N] - 1/l_tild[i]) * prod(lambda[i] - lambda[(i + 1):N])
    }
    logA <- mm * log(M) + sum(log(sapply(1:K, tiny_prod)))
    logrho_k <- 1/(N - K) * sum(log(lambda[(K + 1):N])) - 1/(N - K) * log(sum(lambda[(K +
        1):N]))
    out[K] <- (mm + K)/2 * log(M) - M * (N - K) * (logrho_k) - M/2 * sum(log(lambda[1:K])) +
      M * sum(log(lambda[(K + 1):N]))
    }

     if (evidence) {
      return(out)
    } else {
      return(which.max(out))
    }
}
