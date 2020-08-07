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
