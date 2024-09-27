#' @export
.aValue <- function(fit, top = FALSE, v=5){
  if(top){
    # Calculate a(Psi_0, g), the top value of AEF, where Psi_0 = rho_0'.

    func_1 <- function(w_t){
      res_1 <- (.psi_0(w_t))^2 *.density(w_t, "norm")

      res_1
    }
    func_2 <- function(w_t){
      res_2 <- exp(w_t)/2 *.density(w_t, "norm")

      res_2
    }

    a <- integrate(func_1, -50,50)$value /(integrate(func_2, -50,50)$value)^2

  } else{
    # Calculate a(Psi, g), the bottom value of AEF, where Psi = rho'.
    #topFunc <- function(w){.density(w)*.psi(w)^2}
    #bottomFunc <- function(w){.density.grad(w) * .psi(w)}
    #a <- integrate(topFunc, -50.0, 50.0)$value / (integrate(bottomFunc, -50.0, 50.0)$value)^2

    if(fit$methods == "QML"){
      func_1 <- function(w_t){
        res_1 <- (.psi_0(w_t))^2 *.density(w_t, "norm")

        res_1
      }
      func_2 <- function(w_t){
        res_2 <- exp(w_t)/2 *.density(w_t, "norm")

        res_2
      }

      a <- integrate(func_1, -50,50)$value /(integrate(func_2, -50,50)$value)^2
    } else if(fit$methods == "MLE"){
      func_1 <- function(w_t){
        res_1 <- (.psi(w_t, methods = "QML", dist = "std", v=v))^2 *.density(w_t, "norm")

        res_1
      }
      func_2 <- function(w_t){
        res_2 <-  .psi_p(w_t, methods = "QML", dist = "std", v=v)*.density(w_t, "norm")

        res_2
      }

      a <- integrate(func_1, -50,50)$value /(integrate(func_2, -50,50)$value)^2
    } else{
      func_1 <- function(w_t){
        res_1 <- (.psi(w_t, div_=fit$fixed_pars[1],methods=fit$methods))^2 *.density(w_t, "norm")

        res_1
      }
      func_2 <- function(w_t){
        res_2 <- .psi_p(w_t, div_=fit$fixed_pars[1], methods=fit$methods) *.density(w_t, "norm")

        res_2
      }

      a <- integrate(func_1, -10,10)$value /(integrate(func_2, -10,10)$value)^2
    }
  }

  a
}
.psi <- function(w, div_=.8, methods, dist="norm", v=5){

  if(methods=="QML"){
    p <- .psi_0(w)
  } else if(methods=="MLE"){
    p <- -1/2 + (v+1) /(2 *(1 +(v-2)*exp(-w)))
  } else{

    # psi = rho' = (div * m(rho_0 / div))' = (div * m(-log(g_0(w))/div))'
    # = - m'(-log(g_0(w))/div) * g_0'(w) / g_0(w)

    b <- 4.3 #6.7428
    a <- 4.0 #b-0.5
    ba <- b-a
    x <- (exp(w)-w +log(2*pi))/2 # Note the diff btw x and w.
    x <- x/div_

    u <- as.numeric(x>b)
    v <- as.numeric(x<a)

    p <- v*1 +
      (1-u-v)*(1 -(x-a)^2/ba^2 -(2*a^2*b)/ba^3 +2*(x^3 -(2*a+b)*x^2 +(a^2+2*a*b)*x)/ba^3) +
      u*0

    p <- p * .psi_0(w)
  }

  p
}
# psi'(w_t), derivative of psi(w_t).
.psi_p <- function(w, div_=.8, methods, dist="norm", v=5){

  if(methods=="QML"){
    p_p <- exp(w)/2
  } else if(methods=="MLE"){
    p_p <- ((v+1) *(v-2) *exp(-w)) /(2 *(1 +(v-2)*exp(-w))^2)
  } else{
    # psi' = (m'(-log(g_0(w))/div) *(- g_0'(w) / g_0(w)) )'
    # = m''(-log(g_0(w))/div) *psi_0(w_t)^2 /div
    #   + m'(-log(g_0(w))/div) *psi_0(w_t)'

    b <- 4.3 #6.7428
    a <- 4.0 #b-0.5
    ba <- b-a
    x <- (exp(w)-w +log(2*pi))/2 # Note the diff btw x and w.
    x <- x/div_

    u <- as.numeric(x>b)
    v <- as.numeric(x<a)

    p1 <- v*1 +
      (1-u-v)*(1 -(x-a)^2/ba^2 -(2*a^2*b)/ba^3 +2*(x^3 -(2*a+b)*x^2 +(a^2+2*a*b)*x)/ba^3) +
      u*0
    p2 <- v*0 +
      (1-u-v)*(0 -2*(x-a)/ba^2 -0 +2*(3*x^2 -2*(2*a+b)*x +(a^2+2*a*b))/ba^3) +
      u*0

    p_p <- p1 *exp(w)/2 + p2 *.psi_0(w)^2 /div_
  }

  p_p
}
.psi_0 <- function(w){
  #psi_0 = rho_0' = (-log(g_0))' = - g_0'/g_0
  p <- (exp(w)-1)/2

  p
}
# g(w), the density of w_t.
.density <- function(w, dist="norm", v=5){
  if(dist=="norm"){
    g <- exp((w-exp(w))/2)/sqrt(2*pi)
  } else if(dist=="std"){
    g <- gamma((v+1)/2) *exp(w/2) /(sqrt(v-2) *gamma(v/2) *(1+exp(w)/(v-2))^((v+1)/2))
  } else{
    stop("the current version only apply for the assumption where the density of z_t is symmetric around 0.")
  }

  g
}
# derivative of g(w).
.density.grad <- function(w, dist="norm", v=5){
  if(dist=="norm"){
    g <- exp((w-exp(w))/2) * (1-exp(w))/ (2*sqrt(2*pi))
  } else if(dist=="std"){
    g <- .density(w, "std", v)/2 -.density(w, "std", v) *exp(w) *(v+1) /(2 *(v-2) *(1 +exp(w)/(v-2)))
  } else{
    stop("the current version only apply for the assumption where the density of z_t is symmetric around 0.")
  }

  g
}
#' @export
.expFisherI <- function(fit, true_pars){

  # Big difference between expected and observed Fisher Information matrix.
  T <- length(fit$data)

  if(fit$methods == 'MLE'){
    # to be added.
  } else{

    FS <- matrix(c(0.0, 0.0, 0.0), nrow = 3, ncol = 1, byrow = TRUE)

    FS[1,1] <- 1/(1-true_pars[3])

    # it seems like k is used here as an index
    # could any k was supposed to be the global one implemented before?
    for(k in 1:(T-1)){
      FS[2,1] <- FS[2,1] + true_pars[3]^(k-1) * (fit$data[T-k])^2
    }

    FS[3,1] <- true_pars[1]/(1-true_pars[3])^2
    for(k in 1:(T-2)){
      FS[3,1] <- FS[3,1] + true_pars[2]*k*true_pars[3]^(k-1)*(fit$data[T-k-1])^2
    }

    ht <- true_pars[1]/(1-true_pars[3])
    for(k in 1:(T-1)){
      ht <- ht + true_pars[2] * true_pars[3]^(k-1) * (fit$data[T-k])^2
    }
  }

  FI <- FS %*% t(FS)
  FI <- FI/(2*ht)


  FI
}

