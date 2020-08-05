robGarch11 <- function(x){

  n <- prod(length(x))
  v <- varin(x)
  y <- x/sqrt(v)
  y <- y-median(y)
  for(i in 1:n){
    if(y[i] == 0){
      i <- i
      y[i] <- 10^(-10)
    }
  }
  vini <- varin(y)

  nesb1 <- nnest1(y, vini)
  nesb2 <- nnest2(y, vini)
  mu <- median(x)

  BM1 <- c(nesb1[1]*v, nesb1[2], nesb1[3], mu, KB1)
  BM2 <- c(nesb2[1]*v, nesb2[2], nesb2[3], mu, KB2)

  bms <- c(BM1, BM2)
  return(bms)
}

varin <- function(x){

  n <- prod(length(x))
  a <- tausq(x)
  xc <- x^2
  int <- tausq(xc-1)
  i <- 1
  v <- x[1]^2

  while(v > a+int & i<30){
    i <- i+1
    v <- x[i]^2
  }

  if(i==30){
    error <- 1
  }
  v <- a
  v
}

tausq <- function(x){
  m <- prod(length(x))
  s <- SEST2(x) # Global declaration.
  assign("Sestim", s, envir = .GlobalEnv)
  r <- RHO2(x/s)
  t <- mean(r)*s^2/0.4797
  return(t)
}

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
  return(s)
}

# A continous rho function.
RHO2 <- function(x){

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

  return(ps)
}


nnest1 <- function(y, vini){

  #rm(Muestray)
  #rm(Muestram)
  #rm(Muestrac)

  assign("Muestray", y, envir = .GlobalEnv) # Glocal declaration
  n <- prod(length(y))

  alfa1min <- 0.
  alfa1max <- 1.
  alfa0min <- 0.1
  alfa0max <- 1.
  beta1min <- 0.
  beta1max <- 1.
  nalfa1 <- 5.
  nalfa0 <- 5.
  nbeta1 <- 5.
  ml <- 10^8
  flag <- 0.

  lmalfa1 <- (alfa1max-alfa1min)/nalfa1
  lmalfa0 <- (alfa0max-alfa0min)/nalfa0
  lmbeta1 <- (beta1max-beta1min)/nbeta1
  yc <- y[1:n]^2

  assign("Muestrac", yc, envir = .GlobalEnv) # Glocal declaration
  y2 <- log(yc)
  assign("Muestram", y2, envir = .GlobalEnv) # Glocal declaration
  for(ni in 0:nalfa0){
    for(nj in 0:nalfa1){
      for(nk in 0:nbeta1){

        alfa1 <- alfa1min + nj*lmalfa1
        alfa0 <- alfa0min + ni*lmalfa0
        beta1 <- beta1min + nk*lmbeta1

        var <- rep(0.0, n)
        if(alfa1+beta1<1){
          var[1] <- vini

          for(k in c(5, 20)){ # Question?????
            l <- k+1

            for(i in 2:n){
              var[i] <- alfa0+(alfa1*gk(yc[i-1]/var[i-1],k,l)+beta1)*var[i-1]
            }

            nml <- mean(nfun(y2[2:n]-log(var[2:n])))
            if(nml<ml){
              flag <- 1
              vi <- c(alfa0,alfa1,beta1)
              ml <- nml
            }
          }
        }
      }}}
  # we now optimize
  alfar <- nloptr(c(vi, vini), Fnue) # fminsearch(Fnue,vi,vini)
  nes <- alfar$xmin # c(alfar[1],alfar[1,2],alfar[1,3])

  return(nes)
}

Fnue <- function(vi, vini){
  n <- prod(length(Muestram))
  y2 <- Muestram
  yc <- Muestrac
  var <- rep(0.0, n)
  var[1] <- vini
  nml <- 10^7
  if(vi[1]>0 & vi[2]>= 0 & vi[3]>=0){
    for(k in c(5, 30)){
      l <- k+1
      for(i in 2:n){
        var[i]<-vi[1]+(vi[2]*gk(yc[i-1]/var[i-1],k,l)+vi[3])*var[i-1]
      }
      ml<-mean(nfun(y2[2:n]-log(var[2:n])))
      if(ml<nml){
        nml<-ml
        if(k == 5){
          assign("KB1", 1, envir = .GlobalEnv)
        }
        else{
          assign("KB1", 0, envir = .GlobalEnv)
        }
      }
    }
  }

  return(nml)
}

nfun <- function(x){
  b <- 6.7428
  b1 <- b-0.1
  c <- 1
  x <- exp(x/c)-x/c
  ps<- freg(x, b1, b)
  return(ps)
}

freg <- function(x, a,b){
  n <- prod(length(x))
  u <- as.numeric(x>b)
  v <- as.numeric(x<a)

  c1 <- a-(2/(b-a)^3)*((-1/4)*a^4+(1/3)*(2*a+b)*a^3+(1/2)*(-a^2-2*a*b)*a^2+a^2*b*a)
  c2 <- (-1/(3*(b-a)^2))*(b-a)^3+b-(2/(b-a)^3)*((-1/4)*b^4+(1/3)*(2*a+b)*b^3+(1/2)*(-a^2-2*a*b)*b^2+a^2*b*b)

  g=x*v+(1-u-v)*((-1/(3*(b-a)^2))*(x-a)^3+x-(2/(b-a)^3)*((-1/4)*x^4+(1/3)*(2*a+b)*x^3+(1/2)*(-a^2-2*a*b)*x^2+a^2*b*x)+a-c1)+(c2+a-c1)*u

  return(g)
}

gk <- function(x, k, l){
  n <- prod(length(x))
  bot <- 2*k-2*l
  a <- 1/bot
  b <- -2*l/bot
  c <- k^2/bot

  u <- as.numeric(x>l)
  v <- as.numeric(x<k)

  g<- x*v + (1-u-v)*(a*x^2 + b*x + c) + u*(a*l^2+b*l+c)

  return(g)
}

nnest2 <- function(y, vini){
  rm(Muestray)
  rm(Muestram)
  rm(Muestrac)

  assign("Muestray", y, envir = .GlobalEnv)
  n <- prod(length(y))

  alfa1min <- 0.
  alfa1max <- 1.
  alfa0min <- 0.1
  alfa0max <- 1
  beta1min <- 0.
  beta1max <- 1.
  nalfa1 <- 5
  nalfa0 <- 5
  nbeta1 <- 5
  ml <- 10^8
  flag <- 0
  lmalfa1 <- (alfa1max-alfa1min)/nalfa1
  lmalfa0 <- (alfa0max-alfa0min)/nalfa0
  lmbeta1 <- (beta1max-beta1min)/nbeta1
  yc <- y[1:n]^2

  assign("Muestrac", yc, envir = .GlobalEnv)
  y2 <- log(yc)
  assign("Muestram", y2, envir = .GlobalEnv)

  for(ni in 0:nalfa0){
    for(nj in 0:nalfa1){
      for(nk in 0:nbeta1){
        alfa1 <- alfa1min+nj*lmalfa1
        alfa0 <- alfa0min+ni*lmalfa0
        beta1 <- beta1min+nk*lmbeta1

        var <- rep(0.0, n)
        if(alfa1+beta1 <1){
          var[1] <- vini

          for(k in c(3, 20)){
            l <- k+1
            for(i in 2:n){
              var[i] <- alfa0+(alfa1*gk(yc[i-1]/var[i-1],k,l)+beta1)*var[i-1]
            }

            nml <- mean(nfun2(y2[2:n]-log(var[2:n])))

            if(nml < ml){
              flag <- 1
              vi <- c(alfa0,alfa1,beta1)
              ml <- nml
            }
          }
        }
      }
    }
  }

  alfar <- nloptr(c(vi, vini), Fnue2) #fminsearch(Fnue2,vi,vini)
  nes <- alfar$xmin # c(alfar[1,1],alfar[1,2],alfar[1,3])
}

nfun2 <- function(x){

  b <- 6.7428
  b1 <- b-0.5
  c <- 0.85
  x <- exp(x/c)-x/c
  ps <- freg(x,b1,b)

  return(ps)
}

Fnue2 <- function(vi, vini){
  n <- prod(length(Muestram))
  y2 <- Muestram
  yc <- Muestrac

  var <- rep(0.0, n)
  var[1] <- vini

  nml <- 10^7

  if(vi[1] >0 & vi[2] >=0 & vi[3] >=0){
    for(k in c(3,30)){
      l <- k+1
      for(i in 2:n){
        var[i]<-vi[1]+(vi[2]*gk(yc[i-1]/var[i-1],k,l)+vi[3])*var[i-1]
      }
      ml <- mean(nfun2(y2[2:n]-log(var[2:n])))
      if(ml<nml){
        nml <- ml
        if(k==3){
          assign("KB2", 1, envir = .GlobalEnv)
        }
        else{
          assign("KB2", 0, envir = .GlobalEnv)
        }
      }
    }
  }

  return(nml)
}



