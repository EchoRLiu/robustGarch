#' @title Summary for rg class
#'
#' @description Summary for rg S3 class
#'
#' @param fit A RG fit object of class \code{\link{rg}}
#' @param main_name the title of the plot, default is "Conditional SD (vs returns)"
#' @param original_ a logical argument. If TRUE, the original return will be plotted. Default is FALSE
#' @param pctReturn_ a logical argument. IF TRUE, the plot function will plot the returns in percentage instead of original. Default is TRUE.
#' @param abs_ a logical argument, when TRUE, the plot function will plot abs(returns) with conditional standard deviation instead of returns, default to TRUE.
#' @param dist one of "std_norm" or "t_dist", the distribution of z_t, default is "std_norm"
#'
#' @name rg-summary
#' @aliases summary.rg
#' @aliases print.rg
#' @aliases plot.rg
#' @aliases aef
#'
#' @examples
#'
#' data("rtn")
#' fit <- robGarch(rtn, methods="BM", fixed_pars = c(0.8, 3.0), optimizer="Rsolnp", stdErr_method = "numDeriv")
#' summary(fit)
#' print(fit)
#' plot(fit)
#' aef(fit)
#'
#' @rdname rg-summary
#' @export
summary.rg <- function(fit){

  res <- rbind(fit$fitted_pars, fit$standard_error, fit$t_value, fit$p_value)
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("Fitted Pars", "Standard Error", "T-value", "P-value")

  cat("Model: ", fit$methods, "\n")
  if (fit$methods == "BM"){
    cat("with div = ", fit$fixed_pars[1], ", k = ", fit$fixed_pars[2], "\n")
  }
  if (fit$methods == "M"){
    cat("with div = ", fit$fixed_pars[1], "\n")
  }
  if (fit$methods == "QML" || fit$methods == "MLE"){
    cat("Assumed distribution: ", fit$distribution.model, "\n")
  }
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
  cat("\nLog-likelihood: ", fit$objective)
  cat("\n\nOptimizer: ", fit$optimizer)
  cat("\nStarting point for optimization: ", fit$optimizer_x0)
  cat("\nTime elapsed: ", fit$time_elapsed)
  cat("\nConvergence Message: ", fit$message)
}
#' @rdname rg-summary
#' @export
print.rg <- function(fit){

  res <- rbind(fit$fitted_pars, fit$standard_error, fit$t_value, fit$p_value)
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("Fitted Pars", "Standard Error", "T-value", "P-value")

  cat("Model: ", fit$methods, "\n")
  if (fit$methods == "BM"){
    cat("with div = ", fit$fixed_pars[1], ", k = ", fit$fixed_pars[2], "\n")
  }
  if (fit$methods == "M"){
    cat("with div = ", fit$fixed_pars[1])
  }
  if (fit$methods == "QML" || fit$methods == "MLE"){
    cat("Assumed distribution: ", fit$distribution.model, "\n")
  }
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
}
#' @rdname rg-summary
#' @export
plot.rg <- function(fit, pctReturn_ = TRUE, abs_ = TRUE, original_ = FALSE, main_name = "Conditional Volatility (vs |pctReturns(%)|)"){

  .plot.garchsim(fit, pctReturn_, abs_, original_, main_name)

}
#' @rdname rg-summary
#' @export
aef <- function(fit, v=5){

  aTop <- .aValue(fit, TRUE, v)
  aBottom <- .aValue(fit, FALSE, v)

  aef <- aTop/aBottom

  cat("The AEF of the estimate with respect to QML is: ", aef, "\n")

  aef
}


