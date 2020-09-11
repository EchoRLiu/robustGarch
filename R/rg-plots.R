#' @import xts
#' @import zoo

.plot.garchsim <- function(fit, main_name = "Conditional Volatility (vs |Returns|)", abs_ = TRUE)
{

  n <- length(fit$data) + 100
  a0 <- fit$fitted_pars["alpha_0"]
  a1 <- fit$fitted_pars["alpha_1"]
  b1 <- fit$fitted_pars["beta_1"]
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

  returns <- zoo(fit$data)
  plot.zoo(returns, type = "l", xlab = "", ylab = "Return", main = "Returns")

  sigma2 <- fit$sigma
  if(abs_){
    retAndVol <- cbind(abs(returns),zoo(sqrt(sigma2), order.by = index(returns)))
  } else{
    retAndVol <- cbind(returns,zoo(sqrt(sigma2), order.by = index(returns)))
  }
  plot.zoo(retAndVol, screens = "single", type = "l", xlab = "",ylab = "",main = main_name,lty = c("dotted","solid"),col = c("blue","black"), lwd = c(.8,1.5))
  if(method == "bounded MEst"){
    legend("topleft", bty = "n",
         legend=as.expression(c(bquote(alpha[0]~'='~.(signif(a0, 2))),
                                bquote(alpha[1]~'='~.(signif(a1, 2))),
                                bquote(beta[1]~'='~.(signif(b1, 2))),
                                bquote(div~'='~.(signif(fit$fixed_pars[1], 2))),
                                bquote(c~'='~.(signif(fit$fixed_pars[2], 2))))))
  } else if(method == "modified MEst"){
    legend("topleft", bty = "n",
           legend=as.expression(c(bquote(alpha[0]~'='~.(signif(a0, 2))),
                                  bquote(alpha[1]~'='~.(signif(a1, 2))),
                                  bquote(beta[1]~'='~.(signif(b1, 2))),
                                  bquote(div~'='~.(signif(fit$fixed_pars[1], 2))))))
  } else if(method == "QML"){
    legend("topleft", bty = "n",
           legend=as.expression(c(bquote(alpha[0]~'='~.(signif(a0, 2))),
                                  bquote(alpha[1]~'='~.(signif(a1, 2))),
                                  bquote(beta[1]~'='~.(signif(b1, 2))))))
  } else{
    NA
  }
  legend("topright",legend = c("Returns","CondVol") ,lty = c("dashed","solid"),
         col = c("blue","black"),lwd = c(0.8, 1.5), bty = "n")

}
