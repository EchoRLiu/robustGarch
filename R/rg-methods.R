#' @title Methods for rg class
#'
#' @description Methods fro rg S3 class
#'
#' @param fit A RG fit object of class \code{\link{rg}}
#' @param model one of robust garch "1" or "2" model,
#'
#'
#'
#'
#'
#'
#'
#' @rdname rg-methods
#' @export
summary.rg <- function(fit, std_err = "numDeriv"){

  model <- fit$mdoel
  pars <- rbind(fit$alpha_0, fit$alpha_1, fit$beta_1)

  # TBD.
}
