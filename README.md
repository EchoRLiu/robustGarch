# robustGarch

**robustGarch** is an R package aiming to provide a method for modelling robust Garch processes (RG), addressing the issue of robustness toward additive outliers - instead of innovations outliers. This work is based on Muler and Yohai (2008) (MY).

## Installation

The package can be installed as following:

```r
devtools::install_github("EchoRLiu/robustGarch")
library(robustGarch)
```

## Example

This is a basic example which shows you how to fit your daily return time series data into robust Garch(1,1) model.

```r
data(gspc)
fit <- robGarch(gspc, methods = "BM", fixed_pars = c(0.8, 3.0), optimizer="Rsolnp", stdErr_method = "numDeriv")
summary(fit)
plot(fit)
```

For more examples and explanation, please refer to the  [robustGarch-Vignette](https://github.com/EchoRLiu/robustGarch/blob/master/vignettes/robustGarch_Vignette.pdf).

## Future Development

Any future development will be released in the github page. A few key features will be added to the package in September 2020:
  
 * Fix the issue with singularity error with Hessian matrix
 * Statistics tests such as std_error, t_value, p_value for Garch parameters
 * Code debug on model filter for M model and QML
 * More optimization choices
 * Extension to robust Garch(p, q)
 * Name changes for better collaboration

<!-- badges: start -->
[![R-CMD-check](https://github.com/EchoRLiu/robustGarch/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/EchoRLiu/robustGarch/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
