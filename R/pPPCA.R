#' Penalized Probabilistic PCA
#'
#' The function returns the results of penalized profile
#'    log-likelihood given a matrix of data or a vector of sample
#'    eigenvalues. The data matrix is assumed to follow the decomposition
#'    \eqn{X = WL + \epsilon}, where rows of \eqn{X} are decomposed to a linear projection
#'    in an orthogonal space plus error. The solution finds the
#'    rank of \eqn{W}, which represents some hidden structure in
#'    the data, such that \eqn{X-WL} have independent and
#'    identically distributed components.
#'
#' @param lambda a numerical vector of sample eigenvalues
#' @param Tvotes the number of possible tuning parameter values to be searched
#' @param condition the condition number for the covariance matrix to determine the maximum dimension, the default proportion is set to 1000
#' @param verbose a logical to indicate whether the details of the penalized voting results should be shown
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
#' pPPCA(lambda = lambda) # supply the sample eigenvalues
#' }
#'
#' @author Wei Q. Deng, \email{deng@@utstat.toronto.edu} and Radu V. Craiu \email{craiu@@utstat.toronto.edu}
#'
#' @references Tipping, M. E., and Bishop, C. M. (1999). Probabilistic principal component analysis.
#'   \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#'   \strong{61}(3), 611-622.
#'
#' @keywords probabilistic PCA, penalized probabilistic PCA, profile log-likelihood,
#'    penalty tuning parameter, effective dimension
#'
#' @export
#'

 pPPCA <- function(lambda = lambda , condition = 1000, Tvotes = 5000, verbose = F){


    if (is.null(lambda)) {
        stop("Please provide a numerical vector of sample eigenvalues")
    }

    if (!requireNamespace("Matrix", quietly = TRUE)) {
        stop("Package Matrix is needed for this function to work. Please install it.",
            call. = FALSE)
    }

    if (!requireNamespace("mvtnorm", quietly = TRUE)) {
        stop("Package mvtnorm is needed for this function to work. Please install it.",
            call. = FALSE)
    }
    if (!requireNamespace("pracma", quietly = TRUE)) {
        stop("Package pracma is needed for this function to work. Please install it.",
            call. = FALSE)
    }

 
     n <- sum(lambda > 0)
     lambda[1:n] <- lambda[1:n]*n/sum(lambda[1:n])
     sigma2 <- sapply(1:(n - 1), function(x) sum(lambda[(x + 1):n])/(n - x))
     l_original <- -0.5*(cumsum(log(lambda)[1:(n-1)]) + (n-(1:(n-1)))*log(sigma2) + n);
     N <- max(which(max(abs(lambda), na.rm=T)/abs(lambda) < condition))

delta_bounds1 <- sapply(3:(N-1), function(x) (lambda[x] - sigma2[x])/(x*lambda[x]/(n-x) + sigma2[x])) ## so that lk-l(k-1) and lk-l(k+1) are respectively convex and concave

delta_bounds2 <- n*sapply(3:(N-2), function(x) (1-sigma2[x])*(1/x-1/n))
# so that sigma2P < 1

#delta_bounds3 <- sapply(2:(N-1), function(x) (lambda[x] - sigma2[x])*(n-x)/x/lambda[x])


   delta_min = max(min(delta_bounds1[delta_bounds1>0], na.rm=T), min(delta_bounds2[delta_bounds2 > 0], na.rm=T))
   delta_max = min(max(delta_bounds1, na.rm=T), max(delta_bounds2, na.rm=T))

  delta_choice = sort(c(delta_min, delta_min*exp(log(delta_max/delta_min)/Tvotes)^(1:(Tvotes-1))))

#ppca_log_penal(delta = delta_min, lambda = evals, lk=F)
#ppca_log_penal(delta = delta_max, lambda = evals, lk=F)

  new_result1 <- sapply(1:length(delta_choice), function(d) ppca_log_penal(delta = delta_choice[d], lambda = lambda, lk=F))
  votes <- table(new_result1); votes
  nComp <- as.integer(names(sort(votes, decreasing=TRUE)[1])); nComp

 if (verbose == T) {
 	return(list(delta_choice, new_result1))
 } else{
	return(nComp)
}
}


#'
#' Penalized profile log-likelihood of a PPPCA model
#'
#' The function returns either the penalized profile
#'     log-likelihood or the value that maximizes the penalized
#'     profile log-likelihood for a given tuning parameter value.
#'
#' @param delta the value of the tuning parameter, must be a positive real number
#' @param lambda a numerical vector of positive sample eigenvalues
#' @param lk a logical specifying whether the penalized profile log-likelihood or
#'     the integer that maximizes the penalized profile log-likelihood should be returned.
#' @return an integer K that maximizes the penalized profile log-likelihood for the given \code{delta} value.
#'
#' @author Wei Q. Deng, \email{deng@@utstat.toronto.edu}
#'
#' @keywords probabilistic PCA, penalized probabilistic PCA, profile log-likelihood, penalty
#'
#' @export
#'



 pPPCA_AIC <- function(lambda = lambda , condition = 900, Tvotes = 5000, m = M){

 n <- sum(lambda > 0)
 lambda <- lambda*n/sum(lambda)
 sigma2 <- sapply(1:(n - 1), function(x) sum(lambda[(x + 1):n])/(n - x))
 l_original <- -0.5*(cumsum(log(lambda)[1:(n-1)]) + (n-(1:(n-1)))*log(sigma2) + n);
 N <- max(which(max(abs(lambda), na.rm=T)/abs(lambda) < condition))


delta_bounds1 <- sapply(2:(N-2), function(x) (lambda[x] - sigma2[x])/(x*lambda[x]/(n-x) + sigma2[x])) ## so that lk-l(k-1) and lk-l(k+1) are respectively convex and concave

delta_bounds2 <- n*sapply(2:(N-2), function(x) (1-sigma2[x])*(1/x-1/n))
# so that sigma2P < 1


   delta_min = max(min(delta_bounds1[delta_bounds1>0], na.rm=T), min(delta_bounds2[delta_bounds2 > 0], na.rm=T))
   delta_max = min(max(delta_bounds1, na.rm=T), max(delta_bounds2, na.rm=T))

  delta_choice = sort(c(delta_min, delta_min*exp(log(delta_max/delta_min)/Tvotes)^(1:(Tvotes-1))))

ddAIC <- delta_choice[which.min(sapply(delta_choice, AIC_choose, lambda=lambda, p = m))]; ddAIC
AIC_R <- ppca_log_penal(delta = ddAIC, lambda= lambda, lk=F); AIC_R

return(AIC_R)
}



 pPPCA_BIC <- function(lambda = lambda , condition = 900, Tvotes = 5000, m = M){

 n <- sum(lambda > 0)
 lambda <- lambda*n/sum(lambda)
 sigma2 <- sapply(1:(n - 1), function(x) sum(lambda[(x + 1):n])/(n - x))
 l_original <- -0.5*(cumsum(log(lambda)[1:(n-1)]) + (n-(1:(n-1)))*log(sigma2) + n);
 N <- max(which(max(abs(lambda), na.rm=T)/abs(lambda) < condition))

delta_bounds1 <- sapply(2:(N-2), function(x) (lambda[x] - sigma2[x])/(x*lambda[x]/(n-x) + sigma2[x])) ## so that lk-l(k-1) and lk-l(k+1) are respectively convex and concave

delta_bounds2 <- n*sapply(2:(N-2), function(x) (1-sigma2[x])*(1/x-1/n))
# so that sigma2P < 1


   delta_min = max(min(delta_bounds1[delta_bounds1>0], na.rm=T), min(delta_bounds2[delta_bounds2 > 0], na.rm=T))
   delta_max = min(max(delta_bounds1, na.rm=T), max(delta_bounds2, na.rm=T))

  delta_choice = sort(c(delta_min, delta_min*exp(log(delta_max/delta_min)/Tvotes)^(1:(Tvotes-1))))


ddBIC <- delta_choice[which.min(sapply(delta_choice, BIC_choose, lambda=lambda, p = m))]; ddBIC
BIC_R <- ppca_log_penal(delta = ddBIC, lambda= lambda, lk=F); BIC_R


return(BIC_R)
}



ppca_log_penal <- function(delta, lambda, lk=T){

    lambda <- ifelse(lambda > 0, lambda, 0)
    n <- sum(lambda > 0)
    sigma2 <- sapply(1:(n-1), function(x) sum(lambda[(x+1):n])/(n-x));
    sigma2P <- sapply(1:(n-1), function(x) sigma2[x]*(n-x)/(n-x-x*delta))

    l_original <- -0.5*(cumsum(log(lambda)[1:(n-1)]) + (n-(1:(n-1)))*log(sigma2) + n);

    Kmax <- suppressWarnings(max(which(sigma2P > 0 & sigma2P < 1)))-1
    penal_loglk  <- l_original[1:Kmax] + (- 0.5*((n-(1:Kmax)*(1+delta))*log((n-(1:Kmax))/(n-(1:Kmax)*(1+delta))) - delta*(1:Kmax)*(log(sigma2)[1:Kmax]+1)))

    # penal_loglk <- suppressWarnings(sapply(1:(n-1), function(xx){
     # -0.5*(cumsum(log(lambda))[xx] + (n-xx)*log(sigma2P[xx]) +
      # rev(cumsum(rev(lambda[2:n])))[xx]/sigma2P[xx] + xx) +
      # 0.5*delta*log(sigma2P[xx])*xx
      # }))

    abs_max <- which.max(penal_loglk); #abs_max

    if (Kmax < 3){
    	Kk <- NA
    } else {
    Kk_local <- which(sapply(2:(Kmax-1), function(x) penal_loglk[x] > penal_loglk[x+1] & penal_loglk[x] > penal_loglk[x-1]))+1 ;	Kk_local

    #plot(penal_loglk); abline(v= Kk_local)

    if(length(Kk_local) == 1){
    	Kk <- Kk_local
    } else {
   		Kk <- NA
    }
    }


  if (lk) {
     return(penal_loglk)
     }else{
     return(Kk)
    }
  }

