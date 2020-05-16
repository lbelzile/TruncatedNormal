# TruncatedNormal

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/TruncatedNormal)](https://cran.r-project.org/package=TruncatedNormal)
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-blue.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html) 
[![Downloads](http://cranlogs.r-pkg.org/badges/TruncatedNormal?color=brightgreen)](http://www.r-pkg.org/pkg/TruncatedNormal)

## Truncated multivariate Normal and Student distributions

A collection of functions to deal with the truncated univariate and multivariate normal and Student distributions, described in [Botev (2017)](https://doi.org/10.1111/rssb.12162) and [Botev and L'Ecuyer (2015)](https://doi.org/10.1109/WSC.2015.7408180).

Main features are

- simulation from multivariate truncated Normal and student distributions using an accept-reject algorithm based on minimax exponential tilting.
- (quasi) Monte-Carlo estimation of the distribution function using separation-of-variables together with exponential tilting for provable performances and theoretical upper bound on the error.
- Cholesky decomposition using the reordering algorithm of Gibson, Glasbey and Elston (1994).

To install the latest version from **CRAN**, use 
```R
install.packages("TruncatedNormal")
```
or else install the latest development from Github via

```R
devtools::install_github("lbelzile/TruncatedNormal")
```

