.plot.garchsim <- function(fit, seed = 42, main_name = "Conditional SD (vs returns)", abs_ = TRUE)
{

  n <- length(fit$data) + 100
  a0 <- fit$fitted_pars["alpha_0"]
  a1 <- fit$fitted_pars["alpha_1"]
  b1 <- fit$fitted_pars["beta_1"]

  set.seed(seed)
  z <- rnorm(n)
  x <- rep(0.0, n)
  sigma2 <- rep(0.0, n)
  sigma2[1] <- a0 / (1-b1)
  #x[1] <- fit$data[1]
  #sigma2[1] <- (x[1]/z[1])^2

  for(i in 2:n){
    sigma2[i] <- a0 + a1 * x[i-1]^2 + b1 * sigma2[i-1]
    x[i] <- z[i] * sqrt(sigma2[i])
  }

  sigma2 <- sigma2[101:n]
  par(mfrow = c(2,1), mar=c(2,2,2,2))
  if(abs_){
    plot(abs(fit$data), main=main_name,type="l", col="grey")
  } else{
    plot(fit$data, main=main_name,type="l", col="grey")
    }
  par(new=TRUE)
  lines(sqrt(sigma2), col = "dark blue", main=main_name)
  legend("topleft", bty = "n",
         legend=as.expression(c(bquote(alpha[0]~'='~.(signif(a0, 2))),
                                bquote(alpha[1]~'='~.(signif(a1, 2))),
                                bquote(beta[1]~'='~.(signif(b1, 2))))))

  #x

}
