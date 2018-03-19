# SPAC2
Statistical and Probabilistic Algorithm for Classification and Clustering

This package provides a penalized approach to estimate the number of principal components from either a list of sample eigenvalues or a data matrix on which the sample eigenvalues can be derived in the context of Probabilistic Principal Components Analysis (PPCA). The package also provides various functions to simulate either the sample eigenvalue or the data matrix under specific structures and possibly with violation to normality or independence assumption. The number of PCs uncovered can be then used as inputs for a variety of subsequent analysis such as clustering or classification.

###Quick Start###

The current release is: version 0.9.3 See "release" tab.

To install our R package "SPAC2", you can either run in R directly:

install.packages("devtools") # if you have not installed already devtools::install_github("WeiAkaneDeng/SPAC2")

To load the library in R, simply run in R:

library("SPAC2")
