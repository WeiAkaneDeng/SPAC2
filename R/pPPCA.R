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
#' @param tau  a tolerance threshold for the smallest eigenvalue, the default value is 0.001.
#' @param penalty an integer indicating the type of penalty function to use. The default option is 1, which corresponds to the model in Deng and Craiu (2021).
#' @param beta a numeric between 0 and 1 indicating the weight towards penalty function 1 or 2.
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
#'
#' @keywords probabilistic PCA, penalized probabilistic PCA, profile log-likelihood, penalty tuning parameter, effective dimension
#'
#' @export pPPCA
#'

pPPCA <- function(lambda, Tvotes = 1000, verbose = FALSE, penalty = 1, tau = 0.001, beta = NULL){

  # clean up negative lambda values
  lambda_clean <- ifelse(lambda > tau, lambda, tau);
  n <- length(lambda_clean)
  nn <- sum(lambda_clean > tau)
  sigma2_val <- sapply(1:(nn - 1), function(x) sum(lambda_clean[(x + 1):nn])/(n - x))
  nX <- length(sigma2_val)-1

  if (penalty == 1){

    bounds_ppca1 <- sapply(2:nX, function(ee) (n-ee)*(lambda_clean[ee] - sigma2_val[ee])/((ee)* lambda_clean[ee] + (n-ee)*sigma2_val[ee]))*sqrt(n/nn); #(the G(k) functions in paper for k and k+1)

    delta_min = min(bounds_ppca1, na.rm=T)
    delta_max = max(bounds_ppca1, na.rm=T)

    delta_choice = delta_min*exp(log(delta_max/delta_min)/Tvotes)^(1:(Tvotes))
    new_result1 <- sapply(1:Tvotes, function(d) penal_loglk(delta = delta_choice[d], lam = lambda_clean, penalty=1, n=n))

    votes <- table(new_result1); votes
    nComp <- as.numeric(names(votes)[order(votes, as.numeric(names(votes)), decreasing=TRUE)[1]]); nComp
    #plot(delta_choice, new_result1); abline(h=11)



  } else if (penalty==2) {

    bounds_ppca2 <- sapply(1:nX, function(ee) (lambda_clean[ee+1]-sigma2_val[ee])*(n-ee)/n)
    #(1-sigma2_val)*(n-(1:length(sigma2_val)))/(1:length(sigma2_val))## guarantees that each sigma2t is smaller than 1
    #sapply(1:(length(sigma2_val)-1), function(ee) (lambda[ee+1]-sigma2_val[ee])*(n-ee)/n)
    # G(k) for k and k+1 for 1/sigma2

    delta_min = min(bounds_ppca2, na.rm=T)
    delta_max = max(bounds_ppca2, na.rm=T)

    delta_choice = delta_min*exp(log(delta_max/delta_min)/Tvotes)^(1:(Tvotes))
    new_result1 <- sapply(1:length(delta_choice), function(d) penal_loglk(delta = delta_choice[d], lam = lambda_clean, penalty=2, n = n))

    votes <- table(new_result1); votes
    nComp <- as.numeric(names(votes)[order(votes, as.numeric(names(votes)), decreasing=TRUE)[1]]); nComp



  } else {

    if (is.null(beta)){
      beta0 <- stats::median(sigma2_val)
    }

    bounds_ppca1 <- sapply(1:nX, function(ee) (n-ee)*(lambda_clean[ee] - sigma2_val[ee])/((ee)* lambda_clean[ee] + (n-ee)*sigma2_val[ee]))*sqrt(n/nn); #(the G(k) functions in paper for k and k+1)
    bounds_ppca2 <- sapply(1:nX, function(ee) (lambda_clean[ee+1]-sigma2_val[ee])*(n-ee)/n)
    bounds_ppca3 <- sapply(1:nX, function(ee) min(bounds_ppca1[ee], bounds_ppca2[ee], na.rm=T))

    delta_min = min(bounds_ppca3, na.rm=T)
    delta_max = max(bounds_ppca3, na.rm=T)

    delta_choice = delta_min*exp(log(delta_max/delta_min)/Tvotes)^(1:(Tvotes))
    new_result1 <- sapply(1:length(delta_choice), function(d) penal_loglk(delta = delta_choice[d], lam = lambda_clean, penalty=3, n=n, beta = beta0))

    votes <- table(new_result1); votes
    nComp <- as.numeric(names(votes)[order(votes, as.numeric(names(votes)), decreasing=TRUE)[1]]); nComp

  }


  if (verbose == TRUE) {
    return(list(nComp , delta_choice, new_result1))
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
#' @param tau a tolerance threshold for the smallest eigenvalue, the default value is 0.001.
#' @param M the number of observations used to calculate sample covariance.
#' @param lam a numerical vector of positive sample eigenvalues
#' @param Tvotes the number of possible tuning parameter values to be searched
#' @param AIC a logical indicator to use AIC as the information criterion, if FALSE, then BIC is used; the default option is AIC.
#' @param beta a numeric between 0 and 1 indicating the weight towards penalty function 1 or 2.
#' @return an integer K that maximizes the penalized profile log-likelihood for the given \code{delta} value.
#'
#' @keywords probabilistic PCA, penalized probabilistic PCA, profile log-likelihood, penalty
#'
#' @export pPPCA_INFO
#'

pPPCA_INFO <- function(lam, Tvotes = 100, M, AIC = TRUE, tau = 0.001, beta = NULL){

  lam_clean <- ifelse(lam > tau, lam, tau);
  n <- length(lam_clean)
  nn <- sum(lam_clean > tau)

  lambda <- lam_clean*n/sum(lam_clean)
  sigma2_val <- sapply(1:(nn - 1), function(x) sum(lambda[(x + 1):nn])/(n - x))
  l_original <- -0.5*(cumsum(log(lambda)[1:(nn-1)]) + (n-(1:(nn-1)))*log(sigma2_val) + n);

  bounds_ppca1 <- sapply(1:length(sigma2_val), function(ee) (n-ee)*(lambda[ee] - sigma2_val[ee])/((ee)* lambda[ee] + (n-ee)*sigma2_val[ee]))*sqrt(n/nn); #(the G(k) functions in paper for k and k+1)
  range_bounds <- range(bounds_ppca1, na.rm=T)
  delta_choice = sort(range_bounds[2]*exp(log(range_bounds[1]/range_bounds[2])/Tvotes)^(1:(Tvotes)))


  if (AIC) {
  ddINFO <- delta_choice[which.min(sapply(delta_choice, function(ee) AIC_choose(delta=ee, lambda = lambda, p = M, n = n)))];
  INFO_R <- penal_loglk(delta = ddINFO, lam = lambda, n=n);
  } else {
  ddINFO <- delta_choice[which.min(sapply(delta_choice, function(ee) BIC_choose(delta=ee, lambda = lambda, p = M, n = n)))];
  INFO_R <- penal_loglk(delta = ddINFO, lam = lambda, n=n);

 }
  return(INFO_R)
}




BIC_choose <- function(delta, lambda, p, n){

  K  <- penal_loglk(delta = delta, lam = lambda, penalty = 1, n=n)

  if (is.na(K)){
    return(NA)
  } else {
    -p*loglk(lam = lambda, n= n)[[1]][K] + log(p)*(K*n+(K)/2-(K)^2/2)/2 #K
  }
}


AIC_choose <- function(delta, lambda, n, p){
  K  <- penal_loglk(delta = delta, lam = lambda, penalty = 1, n= n)
  if (is.na(K)){
    return(NA)
  } else {
    -p*loglk(lam = lambda, n= n)[[1]][K] + 2*((K)*n+(K)/2-(K)^2/2)/2 #2*K
  }
}





penal_loglk <- function(delta = 0.1, lam, penalty = 1, n, beta) {

  lam_clean <- ifelse(lam > 1e-5, lam, 1e-5)
  nn <- sum(lam_clean>1e-5)
  sigma2_val <- sapply(1:(nn-1), function(x) sum(lam_clean[(x + 1):nn])/(n - x))

  if (penalty == 1){

    sigma2P = sapply(1:(nn-1), function(xx) sum(lam_clean[(xx+1):nn])/(n-xx-delta*xx))

    if (sigma2P[1] >=1){
      Kmax <- 1
    } else {
      Kmax <- 2
      while (sigma2P[Kmax] >= sigma2P[Kmax+1] & Kmax < nn-1){
        Kmax <- Kmax+1
      }
      ## since there is no certain that the delta input is within the range, we need to
      ## back confirm which maximum K we can search.
    }

    if (Kmax == 1 | sigma2P[1] >=1) {

      penlk_full <- -0.5*(log(lam_clean[1])+sum(lam_clean[2:nn])/sigma2P[1] + 1 + n*log(2))

    } else {
      pen = delta*(1:(Kmax))*log(sigma2P[1:(Kmax)])
      penlk <- sapply(1:(Kmax), function(xx){
        -0.5*(sum(log(lam_clean[1:xx])) + (n-xx)*log(sigma2P[xx]) +
                sum(lam_clean[(xx+1):nn])/sigma2P[xx] + xx + n*log(2))
      })
      penlk_full <- penlk +  0.5*pen; #which.max(penlk_full)

    }

  } else if (penalty == 2){

    sigma2P = sapply(1:(nn-1), function(xx) sum(lam_clean[(xx+1):nn])/(n-xx) + delta*xx/(n-xx));

    if (sigma2P[1] >=1){
      Kmax <- 1
    } else {
      Kmax <- 2
      while (sigma2P[Kmax] >= sigma2P[Kmax+1] & Kmax < nn-1){
        Kmax <- Kmax+1
      }
    }

    if (Kmax == 1 | sigma2P[1] >=1) {

      penlk_full <- -0.5*(log(lam_clean[1])+sum(lam_clean[2:nn])/sigma2P[1] + 1 + n*log(2)) + 0.5*(-delta)/(sigma2P[1])

    } else {

      pen = -delta*(1:(Kmax))/(sigma2P[1:(Kmax)])
      penlk <- sapply(1:(Kmax), function(xx){
        -0.5*(sum(log(lam_clean[1:xx])) + (n-xx)*log(sigma2P[xx]) +
                sum(lam_clean[(xx+1):nn])/sigma2P[xx] + xx + n*log(2))
      })
      penlk_full <- penlk + 0.5*pen; #which.max(penlk_full)
    }

  } else {

    sigma2P = sapply(1:(nn-1), function(xx) (sum(lam_clean[(xx+1):nn])+delta*beta*xx)/(n-xx-xx*delta*(1-beta)));

    if (sigma2P[1] >=1){
      Kmax <- 1
    } else {
      Kmax <- 2
      while (sigma2P[Kmax] >= sigma2P[Kmax+1] & Kmax < nn-1){
        Kmax <- Kmax+1
      }
    }

    if (Kmax == 1 | sigma2P[1] >=1) {
      penlk_full <- -0.5*(log(lam_clean[1])+sum(lam_clean[2:nn])/sigma2P[1] + 1 + n*log(2)) + 0.5*(-delta)/(sigma2P[1])

    } else {

      pen = -delta*beta*(1:(Kmax))/(sigma2P[1:(Kmax)]) + (1-beta)* delta*(1:Kmax)*log(sigma2P[1:(Kmax)])

      penlk <- sapply(1:(Kmax), function(xx){
        -0.5*(sum(log(lam_clean[1:xx])) + (n-xx)*log(sigma2P[xx]) +
                sum(lam_clean[(xx+1):nn])/sigma2P[xx] + xx + n*log(2))
      })
      penlk_full <- penlk + 0.5*pen; #which.max(penlk_full)
    }

  }

  K_est <- which.max(c(-0.5*(n+n*log(2)),penlk_full))
  K_final <- ifelse(K_est==(Kmax+1), NA, K_est-1)
  return(K_final)
}


