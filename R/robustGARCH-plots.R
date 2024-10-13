#' @importFrom graphics legend par
#' @importFrom stats density integrate median nlminb pnorm
.plot.garchsim <- function(fit, digits = 3, estimation_pos = "topleft", line_name_pos = "topright", 
                           par_ = par(no.readonly = TRUE), pctReturn_ = TRUE, abs_ = TRUE, original_ = FALSE, main_name = "Conditional Volatility (vs |pctReturns(%)|)")
{

  n <- length(zoo::coredata(fit$data)) + 100
  a0 <- fit$fitted_pars["gamma"]
  a1 <- fit$fitted_pars["alpha"]
  b1 <- fit$fitted_pars["beta"]
  method <- fit$methods

  #set.seed(seed)
  #z <- rnorm(n)
  #x <- rep(0.0, n)
  #sigma2 <- rep(0.0, n)
  #sigma2[1] <- a0 / (1-b1)
  #x[1] <- fit$data[1]
  #sigma2[1] <- (x[1]/z[1])^2
  #for(i in 2:n){
  #  sigma2[i] <- a0 + a1 * x[i-1]^2 + b1 * sigma2[i-1]
  #  x[i] <- z[i] * sqrt(sigma2[i])
  #}
  #sigma2 <- sigma2[101:n]

  returns <- zoo::coredata(fit$data)
  sigma2 <- xts::xts(sqrt(fit$sigma), order.by = zoo::index(fit$data))

  if(original_){
    zoo::plot.zoo(fit$data, type = "l", xlab = "", ylab = "Return", main = "Returns")
  }

  if(pctReturn_){
    returns <- returns * 100
    sigma2 <- sigma2 * 100
  }

  if(abs_){
    retAndVol <- cbind(abs(returns), sigma2)
  } else {
    retAndVol <- cbind(returns, sigma2)
  }

  # Ensure retAndVol retains the dates from fit$data
  retAndVol <- xts::xts(retAndVol, order.by = zoo::index(fit$data))

  par_
  zoo::plot.zoo(retAndVol, screens = "single", type = "l", xlab = "", ylab = "", main = main_name,
                lty = c("dotted", "solid"), col = c("blue", "black"), lwd = c(.8, 1.5))
  exprss <- c(bquote(alpha[0]~'='~.(signif(a0, digits))),
              bquote(alpha[1]~'='~.(signif(a1, digits))),
              bquote(beta[1]~'='~.(signif(b1, digits))))
  if(fit$methods == "MLE"){
    exprss <- c(exprss, bquote(v~'='~.(signif(fit$fitted_pars["shape"], digits))))
  } #else if(method == "BM"){
    #exprss <- c(exprss,
    #            bquote(div~'='~.(signif(fit$fixed_pars[1], digits))),
    #            bquote(c~'='~.(signif(fit$fixed_pars[2], digits))))
  #} else if(method == "M"){
    #exprss <- c(exprss,
    #            bquote(div~'='~.(signif(fit$fixed_pars[1], digits))))
  #}
  legend(estimation_pos, bty = "n",
         legend=as.expression(exprss))
  if(abs_){
    legend(line_name_pos,legend = c("Absolute Returns","Conditional Vol") ,lty = c("dashed","solid"),
           col = c("blue","black"),lwd = c(0.8, 1.5), bty = "n")
  } else {
    legend(line_name_pos,legend = c("Returns","Conditional Vol") ,lty = c("dashed","solid"),
           col = c("blue","black"),lwd = c(0.8, 1.5), bty = "n")
  }

}
