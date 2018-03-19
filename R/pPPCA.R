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
#' @param x a data matrix with the number of rows to be reduced; only complete columns are used.
#' @param eigenvals a numerical vector of sample eigenvalues
#' @param Tvotes the number of possible tuning parameter values to be searched
#' @param var_tol the tolerance on the empirical error variance to be greater than 0,
#'    the default proportion is set to 0.01
#' @param tol a tolerance level for the EM algorithm to terminate computations.
#' @param maxits the maximum number of iterations in EM to find the principal components
#' @param printComp a logical to indicate whether only the principal components should be printed
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
#' pPPCA(x = t(normdata)) # supply the data matrix
#' pPPCA(eigenvals = eigen_values) # supply the sample eigenvalues
#' }
#'
#' @author Wei Q. Deng, \email{deng@@utstat.toronto.edu}
#'
#' @references Tipping, M. E., and Bishop, C. M. (1999). Probabilistic principal component analysis.
#'   \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#'   \strong{61}(3), 611-622.
#'
#' @keywords probabilistic PCA, penalized probabilistic PCA, profile log-likelihood,
#'    penalty tuning parameter, effective dimension
#'

pPPCA <- function(x = NULL, eigenvals = NULL, Tvotes = 5000,
                  var_tol = 1e-5, tol = 1e-06, maxits = 500,
                  printComp = F, verbose = F) {

    if (is.null(x) && is.null(eigenvals)) {
        stop("Please provide either a data matrix or a numerical vector of sample eigenvalues")
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

    if (is.null(eigenvals)) {

        X <- as.matrix(x[,!apply(x, 2, function(xx) sum(is.na(xx)) > 0)])
        M <- ncol(X)
        S <- stats::cov(scale(t(X)), use = "pairwise.complete")
        SS <- as.matrix(Matrix::nearPD(S)$mat)
        evals <- eigen(SS)$val
        evecs <- eigen(SS)$vectors

        lambda <- ifelse(evals > 0, evals, 0)
        n <- sum(evals > 0)
        sigma2 <- sapply(1:(n - 1), function(x) sum(lambda[(x + 1):n])/(n - x))
        N <- min(ifelse(sum(sigma2 < var_tol) > 0,
                        which(sigma2 < var_tol)[1], n),
                        which(cumsum(lambda)/n > 1 - var_tol)[1],
                        na.rm = T)

    } else {

        lambda <- ifelse(eigenvals > 0, eigenvals, 0)
        n <- sum(eigenvals > 0)
        sigma2 <- sapply(1:(n - 1), function(x) sum(lambda[(x + 1):n])/(n - x))
        N <- min(ifelse(sum(sigma2 < var_tol) > 0,
                        which(sigma2 < var_tol)[1], n),
                        which(cumsum(lambda)/n > 1 - var_tol)[1],
                        na.rm = T)
    }

    ln_penal1 <- function(delta, lk=T){

    		lambda <- ifelse(lambda > 0, lambda, 0)
  	  	n <- sum(lambda > 0)
  	    sigma2 <- sapply(1:(n-1), function(x) sum(lambda[(x+1):n])/(n-x));
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

  	# and if the remaining ones are monotone then take the largest one.
  	# max print is when K exceeds n/(1+delta)
  	  if (lk) {
  	  	 return(penal_loglk)
  	  	 }else{
  	   	 return(Kk)
  	   	}
  	  }

   b = sapply(2:(N-1), function(f) (f-1)*log(sigma2[f-1]) - (f)*log(sigma2[f]));
   c = sapply(2:(N-1), function(d) log(lambda[d]/sigma2[d-1]) + (n-d)*log(sigma2[d]/sigma2[d-1]))

  C = ifelse(c>0, NA, c)
  B = ifelse(b<0, NA, b)

  r1_K <- -C/B;
  r1_K <- ifelse(r1_K > 0, r1_K, 0)

  delta_min = max(min(r1_K, na.rm=T), n*(1/(N-1)-1/N)*(1-sigma2[N-1]))
  delta_max = min(max(r1_K, na.rm=T), (1-1/N)*(1-sigma2[1])*n)

  ## if no delta large enough for lk_max to be positive, we will use the min otherwise use the half of the delta bound by K_max.

  delta_choice = sort(c(delta_min, delta_min*exp(log(delta_max/delta_min)/Tvotes)^(1:(Tvotes-1))))

  ### step 1.2. determine p

  new_result1 <- sapply(1:length(delta_choice), function(d) ln_penal1(delta= delta_choice[d], lk=F))
  # Votes_LK <- tapply(delta_choice, new_result1, range);
  # Del_range <- as.numeric(c(unlist(Votes_LK[length(Votes_LK)-2])[2],  unlist(Votes_LK[2])[1]));
  votes <- table(new_result1)
  #[delta_choice <  Del_range[2] & delta_choice > Del_range[1]]);
  nComp <- as.integer(names(sort(votes, decreasing=TRUE)[1]));

    if (printComp == F) {
        if (verbose == T) {
            if (nComp == 1) {
                warning("Detected only one cluster and no cluster assignments available")
                return(list(nComp = nComp, Votes = new_result1, delta_range = delta_choice))
            } else {
                return(list(nComp = nComp, Votes = new_result1, delta_range = delta_choice))
            }
        } else {
            if (nComp == 1) {
                warning("Detected only one cluster and no cluster assignments available")
                return(list(nComp = nComp, delta_range = c(delta_min, delta_max, Tvotes)))
            } else {
                return(list(nComp = nComp, delta_range = c(delta_min, delta_max, Tvotes)))
            }
        }
    } else {

        #------------------ get results principal components

        U <- evecs[, 1:nComp]
        D <- diag(lambda[1:nComp])
        L <- t(MASS::mvrnorm(M, mu = rep(0, nComp), Sigma = diag(rep(1, nComp))))
        Sigma2 <- sum(evals[(nComp + 1):n])/(n - nComp)
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
            W.new <- SS %*% W.old %*% solve(Psi + solve(MM) %*% t(W.old) %*% SS %*%
                W.old)
            sigma2.new <- 1/n * sum(diag((SS - SS %*% W.old %*% solve(MM) %*% t(W.new))))
            L <- solve(MM) %*% t(W.new) %*% (X)
            loglk1 <- mvtnorm::dmvnorm(t(X - W.new %*% L), mean = matrix(0, n),
                                sigma = diag(sigma2.new, n), log = T)
            loglk2 <- mvtnorm::dmvnorm(t(L), mean = matrix(0, nComp),
                                sigma = diag(1, nComp), log = T)

            loglk <- sum(loglk1) - sum(loglk2)

            it <- it + 1
            converged <- ((max(abs(loglk - loglk.old)) <= tol) &
                        (max(abs(sigma2.new - sigma2.old)) <= tol) &
                        (max(abs(psych::tr(t(W.new) %*% W.new) - psych::tr(t(W.old) %*% W.old))) <= tol))
        }

        W <- pracma::orth(W.new)
        evs <- eigen(stats::cov(t(X) %*% W))
        evecs <- evs$vectors
        W <- W %*% evecs

        if (verbose == T) {
        return(list(W_unbiased = W,
                       nComp = nComp,
                       eigenvals = lambda,
                       votes = new_result1,
                       sigma2 = sigma2.new,
                       delta_range = delta_choice))
        }else{
          return(list(W_unbiased = W,
                    nComp = nComp,
                    eigenvals = lambda,
                    votes = new_result1,
                    sigma2 = sigma2.new))
        }
      }
}
