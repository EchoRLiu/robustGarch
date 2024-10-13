#' @title Summary for robustGARCH class
#'
#' @description Summary for robustGARCH S3 class
#'
#' @param fit A robustGARCH fit object of class \code{\link{robGarch}}
#' @param x Same as fit, for plot.robustGARCH and print.robustGARCH
#' @param object Same as fit, for summary.robustGARCH
#' @param digits the number of digits for print and plot, default is 3.
#' @param main_name the title of the plot, default is "Conditional SD (vs returns)"
#' @param estimation_pos string that determines the legend position that specifies gamma, alpha, beta estimations. Choice of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Default is "topleft".
#' @param line_name_pos string that determines the legend position that specifies the names of lines in the plot. Choice of "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". Default is "topright".
#' @param par_ graphical parameters that can be set, which is in the form of par(...). The default is par(no.readonly = TRUE).
#' @param original_ a logical argument. If TRUE, the original return will be plotted. Default is FALSE
#' @param pctReturn_ a logical argument. IF TRUE, the plot function will plot the returns in percentage instead of original. Default is TRUE.
#' @param abs_ a logical argument, when TRUE, the plot function will plot abs(returns) with conditional standard deviation instead of returns, default to TRUE.
#' @param ... # to be written
#' @param nu # to be written
#' @param v degrees of freedom in a Student's t-distribution.
#'
#' @name robustGARCH-summary
#' @aliases summary.robustGARCH
#' @aliases print.robustGARCH
#' @aliases plot.robustGARCH
#' @aliases aef
#'
#' @examples
#'
#' data("gspc")
#' fit <- robGarch(gspc, methods="BM", fixed_pars = c(0.8, 3.0),
#'                 optimizer="Rsolnp", stdErr_method = "numDeriv")
#' summary(fit)
#' print(fit)
#' plot(fit)
#' aef(fit)
#' 
#' @export
summary.robustGARCH <- function(object, digits = 3, ...){
# summary.robustGARCH <- function(fit, digits = 3){

# notes on why the name change
# Undocumented arguments in Rd file 'robustGARCH-summary.Rd'
#      ‘object’
# Functions with \usage entries need to have the appropriate \alias
#    entries, and all their arguments documented.
#    The \usage entries must correspond to syntactically valid R code.
#    See chapter ‘Writing R documentation files’ in the ‘Writing R
#    Extensions’ manual.
  res <- rbind(round(fit$fitted_pars, digits), round(fit$standard_error, digits), round(fit$t_value, digits), round(fit$p_value, digits))
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("Estimates", "Std. Errors", "t-statistic", "p-value")

  cat("Model: ", fit$methods, " ")
  if (fit$methods == "BM"){
    cat("with div = ", fit$fixed_pars[1], ", k = ", fit$fixed_pars[2])
  }
  if (fit$methods == "M"){
    cat("with div = ", fit$fixed_pars[1])
  }
  cat("\nData: ", fit$data_name, "\n")
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
  cat("\nLog-likelihood: ", fit$objective)
  cat("\n\nOptimizer: ", fit$optimizer)
  if (fit$methods == "MLE"){
    cat("\nInitial parameter estimates: ", fit$optimizer_x0[1:3], fit$optimizer_x0[5])
  } else{
    cat("\nInitial parameter estimates: ", fit$optimizer_x0[1:3])
  }
  cat("\nTime elapsed: ", fit$time_elapsed)
  cat("\nConvergence Message: ", fit$message)
}

#' @rdname robustGARCH-summary
#' @export
print.robustGARCH <- function(x, digits = 3, ...){

  res <- rbind(round(fit$fitted_pars, digits))
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("Estimates (Std. Errors)")

  res[1,1] <- gsub("  ","",paste(res[1,1],'( ', round(fit$standard_error[1], digits), ' )'))
  res[1,2] <- gsub("  ","",paste(res[1,2],'( ', round(fit$standard_error[2], digits), ' )'))
  res[1,3] <- gsub("  ","",paste(res[1,3],'( ', round(fit$standard_error[3], digits), ' )'))
  if (fit$methods == "MLE"){
    res[1,4] <- gsub("  ","",paste(res[1,4],'( ', round(fit$standard_error[4], digits), ' )'))
  }

  cat("Model: ", fit$methods, "\n")
  cat("Data: ", fit$data_name, "\n")
  cat("Result:\n")
  noquote(res)
}

#' @rdname robustGARCH-summary
#' @export
plot.robustGARCH <- function(x, digits = 3, estimation_pos = "topleft", line_name_pos = "topright", par_ = par(no.readonly = TRUE), pctReturn_ = TRUE, abs_ = TRUE, original_ = FALSE, main_name = "Conditional Volatility (vs |pctReturns(%)|)", ...){

  .plot.garchsim(fit, digits, estimation_pos, line_name_pos, par_, pctReturn_, abs_, original_, main_name)

}

#' @rdname robustGARCH-summary
#' @export
aef <- function(fit, nu=5){

  aTop <- .aValue(fit, TRUE, nu)
  aBottom <- .aValue(fit, FALSE, nu)

  aef <- aTop/aBottom

  cat("The AEF of the estimate with respect to QML is: ", aef, "\n")

  aef
}