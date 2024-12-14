#' @title Summary for robustGARCH class
#'
#' @name robustGARCH-summary
#' @aliases summary.robustGARCH
#' @aliases print.robustGARCH
#' @aliases plot.robustGARCH
#' @aliases coef.robustGARCH
#' @aliases aef
#'
#' @description Summary for robustGARCH S3 class
#'
#' @importFrom stats coef
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
#' @param nu degrees of freedom in a Student's t-distribution.
#'
#' @examples
#'
#' data("gspc")
#' fit <- robGarch(gspc, fitMethod="BM", robTunePars = c(0.8, 3.0),
#'                 optChoice="Rsolnp", SEmethod = "numDeriv")
#' summary(fit)
#' print(fit)
#' plot(fit)
#' coef(fit)
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
  fit <- object
  method = fit$fitMethod
  robTunePars = fit$robTunePars
  fitted_pars = fit$fitted_pars

  # model spec
  model_str = paste0("Model: ", method, ", ")
  if (method == "BM") {
    pars_str = paste0("div = ", robTunePars[1], ", k = ", robTunePars[2])
  } else if (method == "M") {
    pars_str = paste0("div = ", robTunePars[1])
  }
  data_str = paste0("Data: ", fit$data_name)
  obs_str = paste0("Observations: ", length(fit$data))
  cat("\n")
  cat(model_str)
  cat(pars_str)
  cat("\n")
  cat(data_str)
  cat("\n")
  cat(obs_str)
  cat("\n")

  # fit result
  significance = sapply(fit$p_value, function(p) {
    if (is.na(p)) ""
    else if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else if (p < 0.1) "."
    else ""
  })
  coefficients <- rbind(round(fitted_pars, digits),
                        round(fit$standard_error, digits),
                        round(fit$t_value, digits),
                        round(fit$p_value, digits))
  colnames(coefficients) <- names(fitted_pars)
  rownames(coefficients) <- c("Estimate", "Std.Error", "t-statistic", "p-value")
  coefficients = t(coefficients)
  res = cbind(coefficients, significance)
  cat("\n")
  cat("Result:")
  cat("\n")
  print.default(format(res, digits = digits), print.gap = 2L, quote = FALSE)
  cat("---")
  cat("\n")
  cat("Signif. codes:  0 '***'' 0.001 '**'' 0.01 '*' 0.05 '.' 0.1 '' 1")
  cat("\n")

  # initial estimate
  if(method == "tMLE") {
    idx_pars = c(1,2,3,5)
  } else {
    idx_pars = 1:3
  }
  initial_estimate = c(fit$initialPars[idx_pars])
  names(initial_estimate) = names(fitted_pars)
  cat("\n")
  cat("Initial parameter estimates:")
  cat("\n")
  print.default(format(initial_estimate, digits = digits), print.gap = 2L, quote = FALSE)

  # optimization
  ll_str = paste0("Log-likelihood: ", format(fit$objective, digits=6))
  opt_str = paste0("Optimizer: ", fit$optChoice)
  time_str = paste0("Time elapsed: ", format(fit$time_elapsed, digits=6))
  conv_str = paste0("Convergence Message: ", fit$message)
  cat("\n")
  cat(ll_str)
  cat("\n")
  cat(opt_str)
  cat("\n")
  cat(time_str)
  cat("\n")
  cat(conv_str)
  cat("\n")

  # summary object
  converged = ifelse(fit$message == 0, TRUE, FALSE) # I'm not sure if this is true
  summary_list = list(
    method = method,
    coefficients = coefficients,
    loglikelihood = fit$objective,
    converged = converged
  )
  invisible(summary_list)
}

#' @rdname robustGARCH-summary
#' @export
print.robustGARCH <- function(x, digits = 3, ...) {

  fit = x
  cat("\n")
  cat("Model: ", fit$fitMethod)
  cat("\n")
  cat("\n")
  cat("Coefficients:")
  cat("\n")
  print.default(format(coef(fit), digits = digits), print.gap = 2L, quote = FALSE)

  invisible(fit)
}

#' @rdname robustGARCH-summary
#' @export
plot.robustGARCH <- function(x, digits = 3, estimation_pos = "topleft", line_name_pos = "topright", par_ = par(no.readonly = TRUE), pctReturn_ = TRUE, abs_ = TRUE, original_ = FALSE, main_name = "Conditional Volatility (vs |pctReturns(%)|)", ...){

  fit <- x
  .plot.garchsim(fit, digits, estimation_pos, line_name_pos, par_, pctReturn_, abs_, original_, main_name)

}

#' @rdname robustGARCH-summary
#' @export
coef.robustGARCH = function(object, ...) {

  object$fitted_pars

}

#' @rdname robustGARCH-summary
aef <- function(fit, nu=5){

  aTop <- .aValue(fit, TRUE, nu)
  aBottom <- .aValue(fit, FALSE, nu)

  aef <- aTop/aBottom

  cat("The AEF of the estimate with respect to QML is: ", aef, "\n")

  aef
}