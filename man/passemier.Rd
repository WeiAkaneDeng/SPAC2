% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/passemier2017.R
\name{passemier}
\alias{passemier}
\title{A bias corrected criterion for selecting number of principal components}
\usage{
passemier(x = NULL, lambda = NULL, M = NULL, constant = 0.05)
}
\arguments{
\item{x}{a data matrix with the number of rows to be reduced; only complete columns are used.}

\item{lambda}{a numeric vector of sample eigenvalues of the covariance matrix of t(\code{x})}

\item{M}{if \code{x} were not supplied, \code{M} should be given as the number of columns of \code{x}.}

\item{constant}{a small prefixed constant and set to the recommended value of 0.05. See Passemier et al., (2017) for details.}
}
\value{
an integer K
}
\description{
The function returns the choice for PCA as a by-product of the bias-corrected residual variance estimate.
}
\examples{
\dontrun{
X <- MASS::mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale(X)))$mat))$val
passemier(lambda = eigen_values, M = 100)
passemier(lambda = eigen_values, M = 5000)
}

}
\references{
Passemier, Damien, Zhaoyuan Li, and Jianfeng Yao (2017). On estimation of the noise variance
	in high dimensional probabilistic principal component analysis. \emph{Journal of the Royal Statistical Society:
	Series B (Statistical Methodology)} \strong{79.1}: 51-67.
}
\keyword{Bayesian}
\keyword{Laplace's}
\keyword{criterion,}
\keyword{evidence}
\keyword{information}
\keyword{log-likelihood,}
\keyword{method,}
\keyword{model}
\keyword{profile}
\keyword{selection,}
