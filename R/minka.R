#' Automatic Choice of Dimensionality for PCA using Bayesian Posterior
#'
#' The function returns the choice dimension for PCA under the PPCA setup using Laplace approximation.
#'
#' @param lambda a numeric vector of sample eigenvalues of the covariance matrix of \code{X} of dimension n by M.
#' @param M  the number of columns of \code{X}.
#' @param tau  a tolerance threshold for the smallest eigenvalue, the default value is 0.001.
#' @param BIC a logical indicating whether the Laplace's method or the BIC approximation should be used.
#' @param verbose a logical specifying whether the posterior evidence or
#'     the integer that minimized the evidence should be returned
#' @return an integer K between 1 and n that maximizes the posterior evidence  by Laplace's method or BIC approximation.
#'
#' @examples
#' \dontrun{
#' X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' minka2001(lambda = eigen_values, M = 100, BIC=TRUE)
#' minka2001(lambda = eigen_values, M = 100, BIC=FALSE)
#' minka2001(lambda = eigen_values, M = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection, Laplace's method, Bayesian evidence
#'
#' @references Minka, T. (2000). Automatic choice of dimensionality for PCA. **Advances in neural information processing systems**, *13*, 598-604. [http://dblp.uni-trier.de/db/conf/nips/nips2000.html#Minka00]
#'
#' @export minka2001
#'

minka2001 <- function(lambda, M, verbose=FALSE, tau=0.001, BIC=FALSE) {

  if (is.null(lambda)) {
    stop("Please provide either a data matrix or a numerical vector of sample eigenvalues")
  }

   if (is.null(M)) {
      stop("Please provide the number of observations or features along with the sample eigenvalues")
    }
  n <- length(lambda > tau)

 	kmax <- min(c(n-1, M-2))
	lambda <- lambda[lambda > ifelse(is.null(tau), 1e-3, tau)]
	n = length(lambda)

	if (!BIC) {

	logediff = NA;
	for (i in 1:(length(lambda)-1)){
	  j = (i+1):length(lambda);
	  logediff[i] = sum(log(lambda[i] - lambda[j])) + (n-length(lambda))*log(lambda[i]);
	}
	cumsum_logediff = cumsum(logediff);

# exact match - safe
#% logediff(i) = sum_{j>i} log(e(i) - e(j))
#logediff = zeros(1,length(e));
#for i = 1:(length(e)-1)
#  j = (i+1):length(e);
#  logediff(i) = sum(log(e(i) - e(j))) + (d-length(e))*log(e(i));
#end
#cumsum_logediff = cumsum(logediff);

	n1 = M-1;
	ehat = (lambda)*M/n1;
	inve = 1/lambda;
	invediff = matrix(rep(inve, length(lambda)), ncol=length(lambda)) - matrix(rep(inve, length(lambda)),nrow=length(lambda), byrow=T)
	invediff[invediff <=0 ] = 1;
	invediff = log(invediff);

	cumsum_invediff = apply(invediff, 2, cumsum);
	row_invediff = t(apply(cumsum_invediff, 1, cumsum));

# exact match
#invediff = repmat(inve,1,length(e)) - repmat(inve',length(e),1);
#invediff(invediff <= 0) = 1;
#invediff = log(invediff);
#% cumsum_invediff(i,j) = sum_{t=(j+1):i} log(inve(t) - inve(j))
#cumsum_invediff = cumsum(invediff,1);
#% row_invediff(i) = sum_{j=1:(i-1)} sum_{t=(j+1):i} log(inve(t) - inve(j))
#row_invediff = cumsum(cumsum_invediff, 2);
#%row_invediff = row_sum(cumsum_invediff);
#% row_invediff(k) = sum_{i=1:(k-1)} sum_{j=(i+1):k} log(inve(j) - inve(i))


	loge = log(ehat);
	cumsum_loge = cumsum(loge);
	cumsum_e = cumsum(ehat);

	dn = length(lambda);
	kmax = length(lambda)-1;
	ks = 1:kmax;
	z = log(2) + (n-ks+1)/2*log(pi) - lgamma((n-ks+1)/2);
	cumsum_z = cumsum(z);

# exact match - safe
#dn = length(e);
#kmax = length(e)-1;
#%dn = d;
#%kmax = min([kmax 15]);
#ks = 1:kmax;
#% the normalizing constant for the prior (from James)
#% sum(z(1:k)) is -log(p(U))
#z = log(2) + (d-ks+1)/2*log(pi) - gammaln((d-ks+1)/2);
#cumsum_z = cumsum(z);


	prb <- NA

	for (i in 1:length(ks)){
	  k = i;
	  v = (cumsum_e[length(cumsum_e)] - cumsum_e[k])/(n-k);
	  prb[i] = -n1/2*cumsum_loge[k] + (-n1*(n-k)/2)*log(v);
	  prb[i] = prb[i] - cumsum_z[k] - k/2*log(n1);
	  #% compute h = logdet(A_Z)
	  h = row_invediff[k] + cumsum_logediff[k];
	  #% lambda_hat(i)=1/v for i>k
	  h = h + (n-k)*sum(log(1/v - inve[1:k]));
	  mm = n*k-k*(k+1)/2;
	  h = h + mm*log(M);
	  prb[i] = prb[i] + (mm+k)/2*log(2*pi) - h/2;
	  #% missing terms added August 21 2008
	  prb[i] = prb[i] + 1.5*k*log(2);
	  prb[i] = prb[i] - 0.5*log(n-k);
	}

}	else {

	    prb <- NA

	    for (i in 1:kmax){

	      k = i
	      e1 = lambda[1:k];
	      e2 = lambda[(k+1):length(lambda)];
	      v = sum(e2)/(n-k)
	      likelihood <- -sum(log(e1))  - (n-k)*log(v)
	      mm = n*k - k*(k+1)/2
	      params = mm + k + n + 1
	      prb[i] = likelihood*M - params*log(M)
	    }

	}


	if (verbose) {
	return(prb)
	} else {
	return(which.max(prb))
	}



}









