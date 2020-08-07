# robGarch

**robGarch** is an R package aiming to provide a method for modelling robust Garch processes (RG), addressing the issue of robustness toward additive outliers - instead of innovations outliers. This work is based on Muler and Yohai (2008) (MY).

## Installation

The package can be installed as following:

```js
devtools::install_github("EchoRLiu/robGarch")
library(robGarch)
```

## Example

This is a basic example which shows you how to fit your daily return time series data into robust Garch(1,1) model.

```js
data(rtn)
bm <- robGarch(rtn, bms=1)
summary(bm)
plot(bm)
```

For more examples and explanation, please refer to the  [robGarch-Vignette](https://github.com/EchoRLiu/robGarch/blob/master/vignettes/robGarch_Vignette.pdf).

## Future Development





