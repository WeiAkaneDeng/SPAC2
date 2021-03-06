% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alternative.R
\name{BIC}
\alias{BIC}
\title{Bayesian Information Criterion for PPCA.}
\usage{
BIC(lambda, M)
}
\arguments{
\item{lambda}{a numeric vector of positive sample eigenvalues that sums to N}

\item{M}{a positive integer for the number of observations or features``}
}
\value{
an integer K between 1 and N that minimizes the BIC.
}
\description{
The function returns the dimension that minimized the BIC based on the
   profile log-likelihood while considering all possible dimensions.
}
\examples{
\dontrun{
library(MASS)
X <- mvrnorm(5000, mu = rep(0,10), Sigma = diag(1,10))
eigen_values <- eigen(as.matrix(Matrix::nearPD(stats::cov(scale((X))))$mat))$val
BIC(lambda = eigen_values, M = 1000)
BIC(lambda = eigen_values, M = 5000)
}
}
\references{
Schwarz, Gideon E. (1978), Estimating the dimension of a model,
   \emph{Annals of Statistics}, \strong{6} (2): 461–464, MR 468014,
   doi:10.1214/aos/1176344136
}
\keyword{criterion,}
\keyword{information}
\keyword{log-likelihood,}
\keyword{model}
\keyword{profile}
\keyword{selection,}
