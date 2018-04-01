#' Bayesian Information Criterion for PPCA.
#'
#' The function returns the dimension that minimized the BIC based on the
#'    profile log-likelihood while considering all possible dimensions.
#'
#' @param lambda a numeric vector of positive sample eigenvalues that sums to N
#' @param M a positive integer for the number of observations or features``
#' @return an integer K between 1 and N that minimizes the BIC.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(5000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale((X))))$mat))$val
#' BIC(lambda = eigen_values, M = 1000)
#' BIC(lambda = eigen_values, M = 5000)
#' }
#' @keywords information criterion, profile log-likelihood, model selection,
#'
#' @references Schwarz, Gideon E. (1978), Estimating the dimension of a model,
#'    \emph{Annals of Statistics}, \strong{6} (2): 461–464, MR 468014,
#'    doi:10.1214/aos/1176344136
#' @export

BIC <- function(lambda, M) {

    lambda <- as.numeric(lambda)
    N <- sum(lambda > 0, na.rm = T)
    lambda <- ifelse(lambda > 0, lambda, 0)
    which.min(-2 * ppcaLog(lambda=lambda, M=M) + (((1:(N - 1)) * N + 1
      - (1:(N - 1)) * ((1:(N - 1)) - 1)/2) * log(M)))
}


#' Akaike Information Criterion for PPCA.
#'
#' The function returns the dimension that minimized the AIC based on the
#'    profile log-likelihood while considering all possible dimensions.
#'
#' @param lambda a numeric vector of positive sample eigenvalues that sums to N
#' @param M a positive integer for the number of observations or features``
#' @return an integer K between 1 and N that minimizes the AIC.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale((X))))$mat))$val
#' AIC(lambda = eigen_values, M = 1000)
#' AIC(lambda = eigen_values, M = 5000)
#' }
#'
#' @keywords information criterion, profile log-likelihood, model selection
#'
#' @references Akaike, H. (1974), A new look at the statistical model identification,
#'    \emph{IEEE Transactions on Automatic Control}, \strong{19} (6): 716–723, MR 0423716,
#'    doi:10.1109/TAC.1974.1100705.
#' @export

AIC <- function(lambda, M) {

    lambda <- as.numeric(lambda)
    N <- sum(lambda > 0, na.rm = T)
    lambda <- ifelse(lambda > 0, lambda, 0)
    which.min((-2 * ppcaLog(lambda=lambda, M=M) + (((1:(N -
        1)) * N + 1 - (1:(N - 1)) * ((1:(N - 1)) - 1)/2)) * 2))
}

#----------------------------------------------------------------------------------


#' Alternative methods based on the detection of elbow in sample eigenvalues.
#'
#' The function returns the dimension that is selected based on various methods that
#'    attempt to detect an ``elbow'' in sample eigenvalues.
#'
#' @param lambda a numeric vector of positive sample eigenvalues that sums to $N$.
#' @param methods a character vector of possible methods to detect an elbow in sample eigenvalues.
#' @return an integer $K$ between 1 and $N$.
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
#' @author Wei Q. Deng, \email{deng@@utstat.toronto.edu}
#'
#' @keywords sample eigenvalues, elbow approach
#' @export
#'
elbowEigen <- function(lambda,
                       methods = c("adjD", "cumD", "varD", "cumlog", "logsigma2")){

  all_methods <- c("adjD", "cumD", "varD", "cumlog", "logsigma2")

  lambda <- as.numeric(lambda)
  N <- sum(lambda > 0, na.rm = T)

  adjD <- which.min(c(lambda[-1]/lambda[-N], 1))
  cumD <- which.max(cumsum(lambda)/(1:N)/((sum(lambda) - cumsum(lambda))/(N - 1:N)))
  varD <- which.max(sapply(1:N, function(x) stats::var(lambda[1:x])))
  cumlog <- which.min(log(cumsum(lambda))/(1:N - cumsum(log(lambda))))
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
#' @param lambda a numeric vector of positive sample eigenvalues that sums to $N$
#' @param M a positive integer for the number of observations or features``
#' @return an integer $K$ between 1 and $N$ that minimizes the LRT $p$-value.
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(t(X))))$mat))$val
#' Lawley.Test(lambda = eigen_values, M = 1000)
#' }
#'
#' @author Wei Q. Deng, \email{deng@@utstat.toronto.edu}
#'
#' @references Lawley, D. N. (1956) Tests of significance for the latent roots of
#'    stats::covariance and correlation matrices. \emph{Biometrika} \strong{43}.1/2: 128-136.
#'
#' @keywords hypothesis testing, likelihood ratio test,
#'
#' @export
#'
Lawley.Test <- function(lambda, M) {

  lambda <- as.numeric(lambda)
  N <- sum(lambda > 0, na.rm = T)
  f = M - 1
  pvalue <- sapply(1:N, function(x) {
    p <- N - x
    logQ = sum(log(lambda[x:N])) - p*log(sum(lambda[x:N]/p))
    constant <- -(f - x - 1/6 * (2*p + 1 + 2/p))
    teststat <- constant * logQ
    stats::pchisq(teststat, p * (p + 1)/2 - 1, lower.tail = F)
  })
  which.min(pvalue)
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
#'
#' @author Wei Q. Deng, \email{deng@@utstat.toronto.edu}
#'
#' @references Tipping, M. E., and Bishop, C. M. (1999). Probabilistic principal component analysis.
#'    \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#'    \strong{61}(3), 611-622.
#'
#' @keywords probabilistic PCA, cross validation, profile log-likelihood, model selection
#'
#' @export
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
