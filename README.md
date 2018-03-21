# SPAC2
Statistical and Probabilistic Algorithm for Classification and Clustering

This package implements a penalized profile log-likelihood criterion to estimate the number of effective dimensions of a data matrix. The data structure is modeled similarly as in Probabilistic Principal Components Analysis (PPCA). The package also provides various functions to simulate either the sample eigenvalues or sample data under specific covariance structures and possibly with violation to normality or independence assumption. The effective dimension or the number of principal components uncovered using our approach can be then used as an input in subsequent analysis such as clustering or classification.

###Quick Start###

The current release is: version 0.9.3 See "release" tab.

To install our R package "SPAC2", you can either run in R directly:

install.packages("devtools") # if you have not installed already

devtools::install_github("WeiAkaneDeng/SPAC2")

To load the library in R, simply run in R:

library("SPAC2")
