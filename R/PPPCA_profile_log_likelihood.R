#' Penalized profile log-likelihood of a PPPCA model
#'
#' The function returns either the penalized profile
#'     log-likelihood or the value that maximizes the penalized
#'     profile log-likelihood for a given tuning parameter value.
#'
#' @param delta the value of the tuning parameter, must be a positive real number
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param lambda a numerical vector of positive sample eigenvalues
#' @param M if \code{x} were not supplied, \code{M} should be given as the number of columns of \code{x}.
#' @param lk a logical specifying whether the penalized profile log-likelihood or
#'     the integer that maximizes the penalized profile log-likelihood should be returned.
#' @return an integer K that maximizes the penalized profile log-likelihood for the given \code{delta} value.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' pppcaProfile(delta = 0, lambda = eigen_values, M = 1000, lk = FALSE)
#' pppcaProfile(delta = 0, lambda = eigen_values, M = 1000, lk = TRUE)
#'}
#' @author Wei Q. Deng, \email{deng@@utstat.toronto.edu}
#'
#' @keywords probabilistic PCA, penalized probabilistic PCA, profile log-likelihood, penalty
#'
#' @export
#'
pppcaProfile <- function(delta, x = NULL, lambda = NULL, M = NULL, lk = FALSE) {

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

  kmax = floor(1/(1/n+delta/n))
K_max <- max(2, min(floor(1/(1/n+delta/n/(1-max(sigma2[kmax], sigma2[n-1], na.rm=T)))), n-1))
# control the maximum K is away from singular solutions of delta
# (and since K_max searches contains K_max + 1)

  l_original <- -0.5*(cumsum(log(lambda)[1:(n-1)]) + (n-(1:(n-1)))*log(sigma2) + n);

  penal_loglk  <- l_original[1:K_max] - 0.5*((n-(1:K_max)*(1+delta))*log((n-(1:K_max))/(n-(1:K_max)*(1+delta))) - delta*(1:K_max)*(log(sigma2)[1:K_max]+1))

  penal_loglk <- ifelse(is.finite(penal_loglk), penal_loglk, NA)

K_choices <- ifelse(which.min(diff(penal_loglk)) > 1, which.min(diff(penal_loglk)) , sum(!is.na(penal_loglk)))

## Jan 19, 2018 - fixed the difference here, if the min difference is greater than 1 then the point before the boundary is selected as the maximum choice.

  # total explained variance by penalized model > 0
  Kk <- ifelse(K_max - 1 > sum(diff(penal_loglk)>0, na.rm=T), which(unlist(sapply(2:K_max, function(xx) penal_loglk[xx] > penal_loglk[xx-1] & penal_loglk[xx] > penal_loglk[xx+1])))+1, which.max(penal_loglk))
  Kk <- min(Kk, max(K_choices), which.max(penal_loglk[1:K_choices]), na.rm=T)

    if (lk) {
        return(penal_loglk)
    } else {
        return(Kk)
    }
}
