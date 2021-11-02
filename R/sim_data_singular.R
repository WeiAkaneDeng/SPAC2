#' Simulate Data from Eigenvalue Structure
#'
#' The function returns either the results of penalized profile
#'    log-likelihood given a matrix of data or a vector of sample
#'    eigenvalues. The data matrix has the following decomposition
#'    \eqn{ X = WL + error}, where the rows of \eqn{X} are linear
#'    projections onto the subspace \eqn{ W } by some arbitrary latent
#'    vector plus error. The solution finds the rank of \eqn{W}, which
#'    represents the hidden structure in the data, such that \eqn{ X-WL }
#'    have independent and identically distributed components.
#'
#' @param N the full dimension of the data
#' @param K the true dimension of the
#' @param M the number of features of observations
#' @param sq_singular a vector of numeric values for the squared singular values.
#'    The other parameters can be skipped if this is supplied along with \eqn{N, K, M}.
#' @param sigma2 a positive numeric between 0 and 1 for the error variance.
#' @param last a positive numeric within reasonable range for the difference
#'    between the Kth eigenvalue and the (\eqn{K}+1)th, a very large difference
#'    might not be possible if \eqn{K} is large.
#' @param trend a character, one of \code{equal}, \code{exponential}, \code{linear} and \code{quadratic}
#'   for the type of trend in squared singular values
#' @param rho a numeric value between 0 and 1 for the amount of auto-correlation, i.e. correlation
#'   between sequential observations or features.
#' @param dist a character specifying the error distribution to be one of \code{norm}, \code{t},
#'    for normal and student's t-distribution, respectively.
#' @param df an integer for the degrees of freedom if \code{dist} == \code{t}
#' @param datamat a logical to indicate whether both the data matrix and
#'    sample eigenvalues or only the sample eigenvalues should be returned
#'
#' @return a list containing the simulated data matrix and
#'    sample eigenvalues or a numerical vector of sample eigenvalues.
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' @importFrom Matrix nearPD
#' @importFrom pracma randortho
#'
#' @examples
#' \dontrun{
#' get_data_singular(N = 200, K = 5, M = 1000, sq_singular = c(5,4,2,1,1))
#' get_data_singular(N = 200, K = 5, M = 1000, sigma2 = 0.2, last= 0.1, trend = "exponential")
#' get_data_singular(N = 200, K = 5, M = 1000, sigma2 = 0.8, last= 0.1, trend = "exponential",
#'    rho = 0.2, df = 5, dist = "t")
#' }
#'
#'
#' @export get_data_singular
#'

get_data_singular <- function(N, K, M, sq_singular = NULL,
                              sigma2 = NULL, last= NULL, trend = NULL,
                              rho = NULL, df = NULL, dist = "norm",
                              datamat = TRUE) {

    if (K <= 0 | M <= 0 | N <= 0){
      stop("Please ensure all of N, K, and M are positive integers")
    }

    if (K >= N | K >= M) {
      stop("Please supply an integer K smaller than both N and M")
    }

  N <- as.integer(N);
  K <- as.integer(K);
  M <- as.integer(M);

  if(is.null(sq_singular)){

    if (sigma2 > 1 | sigma2 <= 0){
      stop("Please supply a sigma2 value between 0 and 1")
    }

    sigma2 <- as.numeric(sigma2)
    d2 <- rep(NA, K)

    if (trend == "equal"){
	      d2[K] <- last
	      remain_var <- N-sigma2*N
	      try(if(remain_var < 0) stop("not enough variance left for the first K-1 eigenvalues"));
	      d2[1:(K-1)] <- remain_var/K

	    } else  if (trend == "linear"){
	      d2[K] <- last
	      remain_var <- N-sigma2*N
	      try(if(remain_var < 0) stop("not enough variance left for the first K-1 eigenvalues"));
	      b = (remain_var - (K-1)*d2[K])/(K*(K-1))*2
	      d2[1:(K-1)] <- d2[K] + b*(K-(1:(K-1)))

	    } else if (trend == "quadratic") {
	      d2[K] <- last
	      remain_var <- N-sigma2*N
	      try(if(remain_var < 0) stop("not enough variance left for the first K-1 eigenvalues"));
	      b = (remain_var - (K-1)*d2[K])/(K*(K-1)*(K-1/2))*3
	      d2[1:(K-1)] <- d2[K] + b*(K-(1:(K-1)))^2

	    } else if (trend == "exponential"){
	      d2[K] <- last
	      remain_var <- N-sigma2*N

	      solve_exp <- function(r){
	        (1-r^(K-1))*d2[K] - r^(K-1)*(1-r)*(remain_var)
	      }
	      r <- stats::uniroot(solve_exp, c(0.0001,1-0.0001))$root
	      d2[2:(K-1)] <- d2[K]/r^(K-2:(K-1));
	      d2[1] <- remain_var - sum(d2[-1])
	      d2 <- sort(d2, decreasing=T)
	    }


  }else{

    if (!is.numeric(sq_singular) | length(sq_singular) != K ){
      stop("Please ensure singular is a numerical vector of length K")
    }

    if (sum(sq_singular) >= N ){
      stop("Please ensure sum of the squared singular values is less than N,
           the total amount of standardized variance")
    }
    sigma2 <- 1 - sum(sq_singular)/N
    d2 = sq_singular
   }

	   K = length(d2)
	   U <- pracma::randortho(N)[,1:K]
	   Lambdaa <- U%*%diag(d2)%*%t(U)
	   LeftTerm <- MASS::mvrnorm(M, mu=rep(0, N), Sigma = Lambdaa, empirical = T) # M by N

   if (is.null(rho)){
    rho = 0
   } else if (rho > 1 | rho < 0){
   stop("rho should be a value between 0 and 1")
   }


       if (dist == "norm"){

	      error <- MASS::mvrnorm(M, mu=rep(0, N),Sigma = diag(sigma2, N))

	      if (rho != 0){
	      errorn_AR <- error
	      errorn_AR[1, ] <- error[1, ]
	      for(m in 2:M){
	        errorn_AR[m, ] <- errorn_AR[(m - 1), ]*rho + error[m, ]
	      }
	      } else {
	      	      errorn_AR <- error
	      }

	      }else{

	        errorT <- sqrt(1/3)*mvtnorm::rmvt(M, sigma = diag(sigma2, N), df = df)
	        errorn_AR <- errorT
	        errorn_AR[1, ] <- errorT[1, ]
	        for(m in 2:M){
	          errorn_AR[m, ] <- errorn_AR[(m - 1), ]*rho + errorT[m, ]
	        }
	      }

	      X <- LeftTerm + errorn_AR
	      #X <- MASS::mvrnorm(M, mu=rep(0, N),Sigma = diag(c(d2+sigma2, rep(sigma2, N-K))))

	      sam_eigen <- eigen(stats::cov(X))$val

    if (datamat == TRUE) {
        return(list(X, sam_eigen))
    } else {
        return(sam_eigen)
    }
}
