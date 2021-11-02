#' Information Criterion for PPCA.
#'
#' The function returns the dimension that minimized the AIC or BIC based on the
#'    profile log-likelihood while considering all possible dimensions.
#'
#' @param lambda a numeric vector of positive sample eigenvalues pf length n
#' @param tau  a tolerance threshold for the smallest eigenvalue, the default value is 0.001.
#' @param AIC a logical indicator to use AIC as the information criterion, if FALSE, then BIC is used; the default option is AIC.
#' @param M a positive integer for the number of observations or features.
#' @param verbose a logical specifying whether the posterior evidence or
#'     the integer that minimized the evidence should be returned
#'
#' @return an integer K between 1 and n that minimizes the AIC or BIC.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale((X))))$mat))$val
#' INFO_ppca(lambda = eigen_values, M = 100)
#' INFO_ppca(lambda = eigen_values, M = 100, AIC = FALSE)
#' INFO_ppca(lambda = eigen_values, M = 5000)
#' INFO_ppca(lambda = eigen_values, M = 5000, AIC = FALSE)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection
#'
#' @references Akaike, H. (1974), A new look at the statistical model identification, **IEEE Transactions on Automatic Control**, *19*(6): 716–723, <doi:10.1109/TAC.1974.1100705>
#'
#' @export INFO_ppca

INFO_ppca <- function(lambda, M, verbose=FALSE, tau = 0.001, AIC = TRUE) {

  n <- length(lambda > tau)

  AIClog <- function(lambda, p, n){
    -p*loglk(lam = lambda, n = n)[[1]] + 2*((1:(n-1))*n+(1:(n-1))/2-(1:(n-1))^2/2)/2
  }

  BIClog <- function(lambda, p, n ){
    -p*loglk(lam = lambda, n = n)[[1]] + log(p)*((1:(n-1))*n+(1:(n-1))/2-(1:(n-1))^2/2)/2

  }

  if (AIC){
    infoLK <- AIClog(lambda,p=M, n = n)
  } else {
    infoLK <- BIClog(lambda,p=M, n = n)
  }
 if (verbose) {
   return(infoLK)
 } else {
   return(which.max(infoLK))
 }
}



#----------------------------------------------------------------------------------


#' Alternative methods based on the detection of elbow in sample eigenvalues.
#'
#' The function returns the dimension that is selected based on various methods that
#'    attempt to detect an ``elbow'' in sample eigenvalues.
#'
#' @param lambda a numeric vector of positive sample eigenvalues of length n
#' @param methods a character vector of possible methods to detect an elbow in sample eigenvalues.
#' @return an integer $K$ between 1 and $n$.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
#' elbowEigen(lambda = eigen_values)
#' elbowEigen(lambda = eigen_values)
#' }
#'
#'
#' @keywords sample eigenvalues, elbow approach
#' @export elbowEigen
#'
elbowEigen <- function(lambda, methods = c("adjD", "cumD", "varD", "cumlog", "logsigma2")){

  all_methods <- c("adjD", "cumD", "varD", "cumlog", "logsigma2")

  lambda <- as.numeric(lambda)
  N <- sum(lambda > 0, na.rm = T)

  adjD <- which.min(c(lambda[-1]/lambda[-N], 1))
  cumD <- which.max(cumsum(lambda)/(1:N)/((sum(lambda)-cumsum(lambda))/(N-1:N)))
  varD <- which.max(sapply(1:N, function(x) stats::var(lambda[1:x])))
  cumlog <- which.min(log(cumsum(lambda))- cumsum(log(lambda)))
  logsigma2 <- which.min(log((N -cumsum(lambda[-N]))/(N-1:(N-1)))*(N-1:(N-1)))

  output <- data.frame(adjD, cumD, varD, cumlog, logsigma2)
  names(output) <- all_methods
  return(output[,all_methods %in% methods])
}


#' Likelihood Ratio Tests for equality of last $n-K$ eigenvalues.
#'
#' The function returns the dimension that minimizes the $p$-value from the
#'   likelihood ratio test (LRT) derived by Lawley (1956).
#'
#'
#' @param lambda a numeric vector of positive sample eigenvalues that sums to $n$
#' @param M a positive integer for the number of observations or features
#' @return an integer $K$ between 1 and $n$ that minimizes the LRT $p$-value.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(c(5,4,2,4,1,1,1,1,1,1)))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(t(X))))$mat))$val
#' LawleyTest(lambda = eigen_values)
#' }
#'
#'
#' @references Lawley, D. N. (1956) Tests of significance for the latent roots of covariance and correlation matrices. **Biometrika** *43*.1/2: 128-136. <doi:10.2307/2333586>
#'
#' @keywords hypothesis testing, likelihood ratio test,
#'
#' @export LawleyTest
#'
LawleyTest <- function(lambda, M) {

   n <- sum(lambda > 0)
  ll_evid <- sapply(1:(n-1), function(x) {
    p <- n - x
    l <- mean(lambda[x:n])
    c <- p-1/6*(2*p+1+2/p)+l^2*sum(1/(lambda[1:x]-l)^2)
    teststat <- c*(p*log(sum(lambda[x:n]/p)) - sum(log(lambda[x:n])))
    stats::pchisq(teststat, p*(p+1)/2-1)
  })

  Lawley <- sum(ll_evid==0) + which.min(ll_evid[ll_evid!=0])

  return(Lawley)
  }


#----------------------------------------------------------------------------------

#' Cross-validation for the PPCA model.
#'
#' The function returns the profile log-likelihood
#'    of the PPCA model at respective MLEs for a specific
#'    choice of $K=k$ evaluated by cross-validation.
#'
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param k an integer indicating the tested dimension, should be between 1 and the number of rows of \code{x}.
#' @param fold an integer indicating the number of folds to be used.
#' @return profile log-likelihood of the remaining fold evaluated at MLEs computed from other folds.
#'
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom mvtnorm dmvnorm
#' @importFrom psych tr
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,50), Sigma = diag(1,50))
#' ppcaCV(k = 5, x = X, fold = 5) # 5-fold cross-validation.
#' }
#'#'
#' @references Tipping, M. E., and Bishop, C. M. (1999). Probabilistic principal component analysis. **Journal of the Royal Statistical Society: Series B (Statistical Methodology)**, *61*(3), 611-622. <doi:10.1111/1467-9868.00196>
#'
#' @keywords probabilistic PCA, cross validation, profile log-likelihood, model selection
#'
#' @export ppcaCV
#'
ppcaCV <- function(k = NULL, x = NULL, fold = 5) {

if (is.null(x)) {
  stop("Please provide a data matrix")
}

if (is.null(k) | is.na(k)) {
  stop("Please specify the dimension to be tested in cross-validation")
}

  X <- as.matrix(x[!apply(x, 2, function(xx) sum(is.na(xx)) > 0),])
  k <- as.integer(k)

  if (k > nrow(X) | k < 1) {
    stop("Please ensure k is an integer between 1 and the number of rows of the data matrix x")
  }

  block_size <- floor(dim(X)[1]/fold)
  logf <- NA
  for (i in 1:fold) {
    test <- X[c(1:block_size + (i - 1) * block_size),]
    train <- X[-c(1:block_size + (i - 1) * block_size),]
    list_ppca <- ppcaMLE(t(train), nComp = k)
    WTrain <- list_ppca$W
    sigma2Train <- list_ppca$sigma2
    logf[i] <- ppcaLog(t(test), param = list(WTrain, sigma2Train), EM = FALSE)
    print(paste("it is the ", i, "th fold,", "testing nComp = ", k, ". Please be patient..."))
  }
  return(logf)
}
