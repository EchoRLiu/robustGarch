#' @title Summary for rg class
#'
#' @description Summary for rg S3 class
#'
#' @param fit A RG fit object of class \code{\link{rg}}
#' @param main_name the title of the plot, default is "Conditional SD (vs returns)"
#' @param abs_ a logical argument, when TRUE, the plot function will plot abs(returns) with conditional standard deviation instead of returns, default to TRUE.
#'
#'
#' @name rg-summary
#' @aliases summary.rg
#' @aliases print.rg
#' @aliases plot.rg
#'
#' @examples
#'
#' data("rtn")
#' fit <- robGarch(rtn, methods="bounded MEst", fixed_pars = c(0.8, 3.0), optimizer="Rsolnp", stdErr_method = "numDeriv")
#' summary(fit)
#' print(fit)
#' plot(fit)
#'
#' @rdname rg-summary
#' @export
summary.rg <- function(fit){

  res <- rbind(fit$fitted_pars, fit$standard_error, fit$t_value, fit$p_value)
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("fitted_pars", "standard_error", "t_value", "p_value")

  cat("Model: ", fit$methods, "\n")
  if (fit$methods == "bounded MEst"){
    cat("with div = ", fit$fixed_pars[1], ", k = ", fit$fixed_pars[2], "\n")
  }
  if (fit$methods == "modified MEst"){
    cat("with div = ", fit$fixed_pars[1], "\n")
  }
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
  cat("\nLog-likelihood: ", fit$objective)
  cat("\n\nOptimizer: ", fit$optimizer)
  cat("\n\n")
  cat("\nTime elapsed: ", fit$time_elapsed)
  cat("\nConvergence Message: ", fit$message)
}

#' @rdname rg-summary
#' @export
print.rg <- function(fit){

  res <- rbind(fit$fitted_pars, fit$standard_error, fit$t_value, fit$p_value)
  colnames(res) <- names(fit$fitted_pars)
  rownames(res) <- c("fitted_pars", "standard_error", "t_value", "p_value")

  cat("Model: ", fit$methods, "\n")
  if (fit$methods == "bounded MEst"){
    cat("with div = ", fit$fixed_pars[1], ", k = ", fit$fixed_pars[2], "\n")
  }
  if (fit$methods == "modified MEst"){
    cat("with div = ", fit$fixed_pars[1])
  }
  cat("Observations: ", length(fit$data), "\n")
  cat("\nResult:\n")
  print(res)
}

#' @rdname rg-summary
#' @export
plot.rg <- function(fit, main_name = "Conditional Volatility (vs |Returns|)", abs_ = TRUE){

  .plot.garchsim(fit, main_name, abs_)

}











