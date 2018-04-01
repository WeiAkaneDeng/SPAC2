#' A bias corrected criterion for selecting number of principal components
#'
#' The function returns the choice for PCA as a by-product of the bias-corrected residual variance estimate.
#'
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param lambda a numeric vector of sample eigenvalues of the covariance matrix of t(\code{x})
#' @param M if \code{x} were not supplied, \code{M} should be given as the number of columns of \code{x}.
#' @param constant a small prefixed constant and set to the recommended value of 0.05. See Bai and Ng (2002s) for details.
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
#' Bai2002(lambda = eigen_values, M = 100)
#' Bai2002(lambda = eigen_values, M = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Bai, Jushan, and Serena Ng (2002). Determining the number of factors in approximate factor models.
#'    \emph{Econometrica} \strong{70.1}: 191-221.
#'
#' @export
#'


Bai2002 <- function(x = NULL, lambda=NULL, M = NULL, constant = 0.05, choice = 1, unbiased = T) {

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

  	N <- length(lambda)

  	g1 <- (N+M)/(N*M)*log(N*M/(N+M))
  	g2 <- (N+M)/(N*M)*log(min(N,M))
  	g3 <- log(min(N,M))/min(N,M)

  	if (choice == 1){
  		gg = g1
  	} else if(choice == 2){
  		gg = g2
  	} else{
  		gg = g3
  	}

    for (K in 1:(N-1)){

  	if (unbiased == T){

  		sigma2 <- sum(lambda[(K+1):N])/(N-K)
  		alpha <- lambda[1:K]
  		cn <- N/(M-1)

  		b_sigma2 <- sqrt(cn/2)*(K + sigma2*sum(1/alpha[1:K]))
  		sigma2_star <- sigma2 + b_sigma2/(N-K)*sigma2*sqrt(2*cn)
  		sigma2_star_m <- sigma2 + b_sigma2*sigma2*sqrt(2*cn)

  	PCsj <- sigma2_star + K*sigma2_star_m*gg

  	}else{

  		 VmF <- sum(lambda[(K+1):N])/(N-K)
  		 # estimate of noise variance if model is PPCA
  	PCsj <- (N*M)^(-1)*VmF + K*VmF*gg
  	}
    out[K] <- PCsj
    }

      return(min(which.min(out), N-1))
}
