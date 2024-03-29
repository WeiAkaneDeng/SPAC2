% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PPCA_MLE_EM.R
\name{ppcaMLE}
\alias{ppcaMLE}
\title{MLEs for the PPPCA model}
\usage{
ppcaMLE(x, nComp = 2, tol = 1e-06, maxits = 100)
}
\arguments{
\item{x}{a data matrix with the number of rows to be reduced; only complete columns are used.}

\item{nComp}{an integer specifying the number of principal components or effective dimension retained.}

\item{tol}{a tolerance level for the EM algorithm to terminate computations.}

\item{maxits}{the maximum number of iterations for EM to converge to the MLEs.}
}
\value{
a list of MLEs: the first item of the list is the loading matrix \eqn{W},
   and the second item of the list is the error variance \eqn{\sigma^{2}} or \code{sigma2}.
}
\description{
The function returns the MLEs for the PPCA model at a given choice of dimension.
 This function is adapted from the one made available by Mark Clark on github
 with the link of the original code under ``See also''.
}
\examples{
\dontrun{
library(MASS)
X <- mvrnorm(1000, mu = rep(0,10), Sigma = diag(1,10))
ppcaMLE(x = t(X), nComp = 5)
}

}
\references{
Tipping, M. E., and Bishop, C. M. (1999). Probabilistic principal component analysis. **Journal of the Royal Statistical Society: Series B (Statistical Methodology)**, *61*(3), 611-622. <doi:10.1111/1467-9868.00196>
}
\seealso{
\url{https://github.com/m-clark/Miscellaneous-R-Code/blob/master/ModelFitting/EM\%20Examples/EM\%20algorithm\%20for\%20ppca.R}
}
\keyword{Estimates}
\keyword{Expectation}
\keyword{Likelihood}
\keyword{Maximization,}
\keyword{Maximum}
\keyword{PCA,}
\keyword{and}
\keyword{probabilistic}
