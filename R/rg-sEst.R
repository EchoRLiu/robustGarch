SEST2 <- function(x){

  b <- 1.625
  emed <- 0.675
  s <- 1 # Starting value
  eps <- 1.
  n <- 1

  m <- median(abs(x))/emed
  x <- x/m
  rho1 <- RHO2(x/0.405)
  a <- mean(rho1)/b
  v <- 1-a # Starting value
  si <- a # Starting value
  rho1 <- RHO2(x/(0.405*si)) # Starting value
  a<- mean(rho1)/b # Starting value
  vi <- 1-a # Starting value
  AUX <- v * vi # Starting value

  while (eps> 0.005 & AUX > 0){
    n <- n+1
    s <- si
    v <- vi
    si <- a*s
    rho1 <- RHO2(x/(0.405 * si))
    a <- mean(rho1)/b
    vi <- 1-a
    AUX <- v*vi
    eps <- abs(s-si)/s
  }

  nsec <- 0
  while(eps>0.005){
    ns <- (s+si)/2
    rho1 <- RHO2(x/(0.405*ns))
    a <- mean(rho1)/b
    nv <- 1-a
    AUX <- nv * vi
    if(AUX<0){
      v <- nv
      s <- ns
    }
    else{
      vi <- nv
      si <- ns
    }
    eps<- abs(s-si)/s
    nsec <- nsec + 1
    n <- n+1
  }
  s <- s*m
  if(n>30){
    n <-n
  }

  s
}

# A continous rho function.
RHO2 <- function(x)
{

  n <- prod(length(x))
  G1 <- -1.944
  G2 <- 1.728
  G3 <- -0.312
  G4 <- 0.016

  ax <- abs(x)
  u <- as.numeric(ax>3.0)
  v <- as.numeric(ax<2.0)
  w <- (1-u)*(1-v)
  ps <- v*x^2/2 + w*(G4*x^8/8 + G3*x^6/6 + G2*x^4/4 + G1*x^2/2 + 1.792) + 3.25*u

  ps
}
