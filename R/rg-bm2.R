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

  nes
}

nfun2 <- function(x){

  b <- 6.7428
  b1 <- b-0.5
  c <- 0.85
  x <- exp(x/c)-x/c
  ps <- freg(x,b1,b)

  ps
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

  nml
}
