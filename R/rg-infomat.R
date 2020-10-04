#' @import stats
#' @import numDeriv
#' @export
.aValue <- function(fit, dist = "std_norm", top = FALSE){
  if(top){
    # Calculate a(Psi_0, g), the top value of AEF, where Psi_0 = rho_0'.
    bottomFunc <- function(w){.density.grad(w)*.psi_0(w)}
    a <- -1/integrate(bottomFunc, -50.0, 50.0)$value
  } else{
    # Calculate a(Psi, g), the bottom value of AEF, where Psi = rho'.
    topFunc <- function(w){.density(w)*.psi(w)^2}
    bottomFunc <- function(w){.density.grad(w) * .psi(w)}
    a <- integrate(topFunc, -50.0, 50.0)$value / (integrate(bottomFunc, -50.0, 50.0)$value)^2
  }

  a
}
.psi <- function(w){

  if(method == "QML" || a == b){

    # rho_0'.
    p <- .psi_0(w)
  } else{
    # psi = rho' = (div * m(rho_0 / div))' = (div * m(-log(g_0(w))/div))'
    # = - m'(-log(g_0(w))/div) * g_0'(w) / g_0(w)

    b <- 6.7428
    a <- b-0.5
    ba <- b-a
    x <- (exp(w)-w +log(2*pi))/2 # Note the diff btw x and w.
    x <- x/div

    u <- as.numeric(x>b)
    v <- as.numeric(x<a)

    p <- v*1 +
      (1-u-v)*(1 -(x-a)^2/ba^2 -(2*a^2*b)/ba^3 +2*(x^3 -(2*a+b)*x^2 +(a^2+2*a*b)*x)/ba^3) +
      u*0

    p <- p * .psi_0(w)
    }

  p
}
.psi_0 <- function(w){
  #psi_0 = rho_0' = (-log(g_0))' = - g_0'/g_0
  p <- (exp(w)-1)/2

  p
}
# g(w), the density of w_t.
.density <- function(w, dist="std_norm"){
  if(dist=="std_norm"){
    g <- exp((w-exp(w))/2)/sqrt(2*pi)
  } else if(dist=="t_dist"){
    # add where f is not special norm case, but according to the data.
    stop("under development.")
  } else{
    stop("the current version only apply for the assumption where the density of z_t is symmetric around 0.")
  }

  g
}
# derivative of g(w).
.density.grad <- function(w, dist="std_norm"){
  if(dist=="std_norm"){
    g <- exp((w-exp(w))/2) * (1-exp(w))/ (2*sqrt(2*pi))
  } else if(dist=="t_dist"){
    # add where f is not special norm case, but according to the data.
    stop("under development.")
  } else{
    stop("the current version only apply for the assumption where the density of z_t is symmetric around 0.")
  }

  g
}

