
#-------------------------------------------------------------------
#' MLEs for the PPPCA model
#'
#'  The function returns the MLEs for the PPCA model at a given choice of dimension.
#'  This function is adapted from the one made available by Mark Clark on github
#'  with the link of the original code under ``See also''.
#'
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param nComp an integer specifying the number of principal components or effective dimension retained.
#' @param tol a tolerance level for the EM algorithm to terminate computations.
#' @param maxits the maximum number of iterations for EM to converge to the MLEs.
#'
#' @return a list of MLEs: the first item of the list is the loading matrix \eqn{W},
#'    and the second item of the list is the error variance \eqn{\sigma^{2}} or \code{sigma2}.
#'
#' @importFrom Matrix nearPD
#' @importFrom psych tr
#'
#' @examples
#' \dontrun{
#' library(MASS)
#' X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
#' ppcaMLE(x = t(X), nComp = 5)
#' }
#'
#' @references Tipping, M. E., and Bishop, C. M. (1999). Probabilistic principal component analysis.
#'    \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#'    \strong{61}(3), 611-622.
#'
#' @seealso \url{https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM\%20Examples/EM\%20algorithm\%20for\%20ppca.R}
#'
#' @keywords probabilistic PCA, Expectation and Maximization, Maximum Likelihood Estimates
#'


ppcaMLE <- function(x, nComp = 2, tol = 1e-06, maxits = 100) {


if (is.null(x)) {
  stop("Please provide a data matrix")
}



  X <- x[, !apply(x, 2, function(xx) sum(is.na(xx)) > 0)]
  M <- ncol(X)
  n <- nrow(X)

  if (M < n){
  stop("Please make sure the number of columns exceed the number of rows.")
  }

  S <- stats::cov(scale(t(X)), use = "pairwise.complete")
  SS <- as.matrix(Matrix::nearPD(S)$mat)
  evals <- eigen(SS)$val
  evecs <- eigen(SS)$vectors
  lambda <- evals[evals > 0]
  N <- length(lambda)

  if (nComp == 1) {
    return(list(W = matrix(0, N,N), sigma2 = 1))
  } else{
    U <- evecs[, 1:nComp]
    D <- diag(lambda[1:nComp])
    L <- t(MASS::mvrnorm(M, mu = rep(0, nComp), Sigma = diag(rep(1, nComp))))
    Sigma2 <- sum(lambda[(nComp + 1):N])/(N - nComp)
    W <- U %*% chol(D - Sigma2 * diag(nComp))

    it <- 0
    converged <- FALSE
    loglk <- 0

    while ((!converged) & (it < maxits)) {
      if (exists("W.new")) {
        W.old <- W.new
        sigma2.old <- sigma2.new
      } else {
        W.old <- W
        sigma2.old <- Sigma2
      }

      loglk.old <- loglk
      Psi <- sigma2.old * diag(nComp)
      MM <- t(W.old) %*% W.old + Psi
      W.new <- SS %*% W.old %*% solve(Psi + solve(MM) %*% t(W.old) %*% SS %*% W.old)
      sigma2.new <- 1/N * sum(diag((SS - SS %*% W.old %*% solve(MM) %*% t(W.new))))
      L <- solve(MM) %*% t(W.new) %*% (X)
      loglk1 <- mvtnorm::dmvnorm(t(X - W.new %*% L),
                                 mean = matrix(0, N),
                                 sigma = diag(sigma2.new, N), log = T)
      loglk2 <- mvtnorm::dmvnorm(t(L), mean = matrix(0, nComp),
                                 sigma = diag(1, nComp), log = T)
      loglk <- -sum(loglk1) - sum(loglk2)

      it <- it + 1
      converged <- ((max(abs(loglk - loglk.old)) <= tol) &
                      (max(abs(sigma2.new - sigma2.old)) <= tol) &
                      (max(abs(psych::tr(t(W.new) %*% W.new) - psych::tr(t(W.old) %*% W.old))) <= tol))
    }

    return(list(W = W.new, sigma2 = sigma2.new))
  }
}
