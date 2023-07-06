## frailtyMMpen: Package for Penalized Frailty Models

[![R-CMD-check](https://github.com/heilokchow/frailtyMMpen/workflows/R-CMD-check/badge.svg)](https://github.com/heilokchow/frailtyMMpen/actions)
[![CRAN](http://www.r-pkg.org/badges/version/frailtyMMpen)](https://cran.r-project.org/package=frailtyMMpen)
[![Downloads](https://cranlogs.r-pkg.org/badges/frailtyMMpen)](https://cran.r-project.org/package=frailtyMMpen)
[![Last Commit](https://img.shields.io/github/last-commit/heilokchow/frailtyMMpen)](https://github.com/heilokchow/frailtyMMpen)

This package implements the MM algorithm for a variety types of frailty models which can handle clustered data, multi-event data and recurrent data in addition to the simple frailty model. Besides, this package can obtain the estimation of parameters for penalized regression using LASSO, MCP and SCAD penalties. Currently supported frailty distributions include gamma, log-normal, inverse gaussian and PVF (1<p<2). The estimation procedure is computationally efficient which makes it also capable for handling high-dimensional data. 


Installation
------------

You can install developed version of frailtyMMpen from github with:

``` r
# install.packages("devtools")
devtools::install_github("heilokchow/frailtyMMpen")
```

Example
-------

This is a basic example which shows you how to use this package, you may refer to the package manual for detailed descriptions and examples for each function.

We use the simulated data with 50 clusters and 10 objects in each cluster: 

``` r
data(simdataCL)
```

We first run the non-penalized regression with Gamma frailty and obtain the summary statistics and the plot of conditional baseline hazard.

```r
gam_cl = frailtyMM(Surv(time, status) ~ . + cluster(id), simdataCL, frailty = "gamma")

summary(gam_cl)

plot(gam_cl)

```

Then, we perform the penalized regression with Gamma frailty and LASSO penalty and obtain BIC, degree of freedom under a sequence of tuning parameters and the plot of regularization path.

```r
gam_cl_pen = frailtyMMpen(Surv(time, status) ~ . + cluster(id), simdataCL, frailty = "gamma")

print(gam_cl_pen)

plot(gam_cl_pen)

```

