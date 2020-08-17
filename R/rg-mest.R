#' @title Robust Estimates for GARCH(1,1) Model
#'
#' @name rgMEst
#' @aliases rgmest
#'
#' @description Methods for fitting a Garch(1,1) model with daily log return time series, using two classes of robust extended M-Estimates:
#' (1) maximizing modified likelihood function;
#' (2) M-Estimates with bounds to the propagation of the effect of one outlier on the subsequent predictors of the conditional variances.
#' returns a \code{rg} object
#'
#' @param data a time series of log returns
#' @param class robust M-Estimate class used for Garch(1,1) model, one of "modified MEst", "bounded MEst", default is "bounded MEst"
#' @param class_pars a named numeric vector of parameters of the class. For the "modified MEst" class, the parameter should be "div", which is used in m(rho0(x)) = div * m(rho0(x)/div) and defined as 0 <= div <= 1.0, default is div = 0.85; for the "bounded MEst" class, the parameter should be "div" and "k", which also influences the robustness and defined as k >= 0, default is div = 0.85, k = 3. The smaller div or k is, the more robust the moddel is.
#' @param optimizer optimizer used for optimization, one of "fminsearch", "nloptr", "Rsolnp", default is "fminsearch"
#' @param optimizer_control list of control arguments passed to the optimizer
#' @param stdErr_method method used to calculate standard error, one of "numDerive", "optim", "sandwich", default is "numDeriv" using hessian from numDeriv
#'
#' @return
#' A \code{rg} object(S3), the components of the object are:
#'     \item{data}{the log returns data object for the rg model to be fitted}
#'     \item{classes}{the class called}
#'     \item{class_pars}{named numeric vector of parameters for the class}
#'     \item{optimizer}{the optimizer called}
#'     \item{optimizer_control}{the list of control arguments passed to the optimizer}
#'     \item{optimizer_result}{output of the called optimizer}
#'     \item{stdErr_method}{the method called to calculate standard error}
#'     \item{fitted_pars}{Garch(1,1) parameter estimations output of the called class}
#'     \item{objective}{the optimal likihood value obtained by the optimizer}
#'     \item{time_elapsed}{the time used for the optimization routine}
#'     \item{message}{the message of the convergence status produced by the called solver}
#'     \item{standard_error}{standard erros of the fitted parameters using the method called}
#'     \item{t_value}{t-values of alpha_0, alpha_1, beta_1}
#'     \item{p_value}{p-values of alpha_0, alpha_1, beta_1}
#'
#' @details
#' The \code{rgMEst} function fits a Garch(1, 1) model to a time series of log return data, using one of the two classes of robust extended M-Estimates with certain parameters specified by the user, with guidance and examples from the vignette. The user can also specify the optimizer used during optimization procesure, and the method used to calculate standard error for the fitted parameters.
#'
#' For details of the list of control arguments, please refer to \code{pracma::fminsearch}, \code{nloptr::nloptr}, \code{Rsolnp::solnp}
#'
#' @references Muler, Nora & Yohai, Victor. (2008). Robust estimates for GARCH models. Journal of Statistical Planning and Inference. 138. 2918-2940.
#'
#' @examples
#' # TBD.
#'
#'
#'
#'
#'
#'
#'
#'
#' @rdname rgMEst
#' @export
# Garch(1,1) model fit function
rgMEst <- function(data, classes = c("bounded MEst", "modified MEst" ), class_pars = c(0.85, 3) , optimizer = c("fminsearch", "nloptr", "Rsolnp"), optimizer_control = list(), stdErr_method = c("numDeriv", "optim", "sandwich")){

  if(!is.numeric(data) || length(data)==0)
    stop("Data must be a numeric vector of non-zero length")

  classes = match.arg(classes)
  optimizer = match.arg(optimizer)
  stdErr_method = match.arg(stdErr_method)
  fit <- rgFit_local(data, classes, class_pars, optimizer, optimizer_control)
  # std_err calculation
  if(stdErr_method == "numDeriv")
  {

    H <- numDeriv::hessian(Fnue, x = fit$optimizer_result$pars)/length(data)
    standard_error <- sqrt(diag(abs(solve(H)/length(data))))[1:3]
    t_value <- fit$fitted_pars / standard_error
    p_value <- 2*(1-pnorm(abs(t_value)))

  }
  if(stdErr_method == "optim")
  {
    H <- optim(par = fit$optimizer_result$pars, fn = Fnue,
               method="L-BFGS-B",
               lower=c(-1, -1, -1, 0, class_pars),
               upper=c(1, 1, 1, 2, class_pars+0.00001),
               hessian=TRUE)$hessian/length(data)
    standard_error <- sqrt(diag(abs(solve(H)/length(data))))[1:3]
    t_value <- fit$fitted_pars / standard_error
    p_value <- 2*(1-pnorm(abs(t_value)))

  }
  if(stdErr_method == "sandwich")
  {
    Scores <- matrix(numDeriv::grad(Fnue, x = fit$optimizer_result$pars), nrow = 1)
    V <- (t(Scores) %*% Scores)/length(data)
    H <- optim(par = fit$optimizer_result$pars, fn = Fnue,
               method="L-BFGS-B",
               hessian=TRUE,
               control=list(maxit=1))$hessian/length(data)
    S <- solve(H) %*% V %*% solve(H)

    standard_error <- sqrt(diag(abs(S/length(data))))[1:3]
    t_value <- fit$fitted_pars / standard_error
    p_value <- 2*(1-pnorm(abs(t_value)))
  }

  names(standard_error) <- c("alpha_0", "alpha_1", "beta_1")
  names(t_value) <- c("alpha_0", "alpha_1", "beta_1")
  names(p_value) <- c("alpha_0", "alpha_1", "beta_1")

  fit$standard_error <- standard_error
  fit$t_value <- t_value
  fit$p_value <- p_value

  structure(fit, class="rg")
}


rgFit_local <- function(data, classes, class_pars, optimizer, optimizer_control){

  start_time <- Sys.time()
  # optimizer/optimizer_control

  # Robust Modified M-Estimates Class.
  if(classes == "modified MEst"){
    div <- class_pars[1]
    stop("Modified MEst usage is under development, please use bounded MEst instead.")
  }
  # Robust Bounded M-Estimates Class.
  if(classes == "bounded MEst"){

    n <- length(data)
    v_data <- varin(data)
    data_normalized <- data/sqrt(v_data)
    data_normalized <- data_normalized-median(data_normalized)
    for(i in 1:n){
      if(data_normalized[i] == 0){
        data_normalized[i] <- 10^(-10)
      }
    }
    vini <- varin(data_normalized)

    div <- class_pars[1]
    k <- class_pars[2]
    res <- nnest(data_normalized, vini, div, k, optimizer, optimizer_control)
    optimizer_result <- res
    if(optimizer == "fminsearch"){
      stop("This feature is under active development, please use Rsolnp.")
      fitted_pars <- res$fmin[1:3]
      fitted_pars[1] <- fitted_pars[1]*v_data
      names(fitted_pars) <- c("alpha_0", "alpha_1", "beta_1")
      objective <- res$xopt
      message <- NULL
    } else if(optimizer == "Rsolnp" || optimizer == "nloptr"){
      fitted_pars <- res$pars[1:3]
      fitted_pars[1] <- fitted_pars[1]*v_data
      names(fitted_pars) <- c("alpha_0", "alpha_1", "beta_1")
      objective <- res$objective
      message <- res$message
    } else{
      NA
    }

  }
  time_elapsed <- Sys.time() - start_time

  list(data=data,
       classes = classes,
       class_pars=class_pars,
       optimizer=optimizer,
       optimizer_control=optimizer_control,
       optimizer_result=optimizer_result,
       fitted_pars = fitted_pars,
       objective=objective,
       time_elapsed=time_elapsed,
       message=message)
}

nnest <- function(y, vini, div, k, optimizer, optimizer_control){

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

          for(ki in c(k, 20)){
            l <- ki+1
            for(i in 2:n){
              var[i] <- alfa0+(alfa1*gk(yc[i-1]/var[i-1],ki,l)+beta1)*var[i-1]
            }

            nml <- mean(nfun(y2[2:n]-log(var[2:n]), div))

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

  lb = c( -1.0, -1.0, -1.0, 0.0, div, k)
  ub = c( 1.0, 1.0, 1.0, 1.0, div, k)
  # Limits for optimization, needs to be added. TBD.
  if (optimizer == "fminsearch"){

    res <- pracma::fminsearch(Fnue,c(vi,vini, div, k),lower=lb,upper=ub, method="Hooke-Jeeves")
    return(res)

  } else if (optimizer == "nloptr"){

    stop("This feature is under development, please use Rsolnp instead.")
    res <- nloptr::nloptr(x0 = c(vi, vini, div, k),
                          eval_f = Fnue,
                          lb = lb,
                          ub = ub,
                          opts=optimizer_control)

    return(res)

  } else if (optimizer == "Rsolnp"){

    res <- Rsolnp::solnp(pars = c(vi, vini, div, k),
                         fun = Fnue,
                         LB = lb,
                         UB = ub,
                         control = optimizer_control)

    return(res)

  } else{
    NA
  }
}

varin <- function(x){

  n <- length(x)
  tausq_x <- tausq(x)
  x_square <- x^2
  tausq_xsquare <- tausq(x_square-1)

  i <- 1
  v <- x[1]^2

  while(v > tausq_x+tausq_xsquare & i<30){
    i <- i+1
    v <- x[i]^2
  }

  if(i==30){
    error <- 1
  }
  v <- tausq_x

  v
}

tausq <- function(x){

  s <- sest(x)
  # Global declaration.
  assign("Sestim", s, envir = .GlobalEnv)
  r <- rho(x/s)
  t <- mean(r)*s^2/0.4797

  t
}


sest <- function(x){

  b <- 1.625
  emed <- 0.675
  s <- 1 # Starting value
  eps <- 1.
  n <- 1

  m <- median(abs(x))/emed
  x <- x/m
  rho1 <- rho(x/0.405)
  a <- mean(rho1)/b
  v <- 1-a # Starting value
  si <- a # Starting value
  rho1 <- rho(x/(0.405*si)) # Starting value
  a<- mean(rho1)/b # Starting value
  vi <- 1-a # Starting value
  AUX <- v * vi # Starting value

  while (eps> 0.005 & AUX > 0){
    n <- n+1
    s <- si
    v <- vi
    si <- a*s
    rho1 <- rho(x/(0.405 * si))
    a <- mean(rho1)/b
    vi <- 1-a
    AUX <- v*vi
    eps <- abs(s-si)/s
  }

  nsec <- 0
  while(eps>0.005){
    ns <- (s+si)/2
    rho1 <- rho(x/(0.405*ns))
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
rho <- function(x){

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

nfun <- function(x, div){

  b <- 6.7428
  b1 <- b-0.5
  x <- exp(x/div)-x/div
  ps <- freg(x,b1,b)

  ps
}

Fnue <- function(start_pars){

  n <- prod(length(Muestram))
  y2 <- Muestram
  yc <- Muestrac

  vi <- start_pars[1:3]
  vini <- start_pars[4]
  div <- start_pars[5]
  k <- start_pars[6]

  var <- rep(0.0, n)
  var[1] <- vini

  nml <- 10^7

  if(vi[1] >0 & vi[2] >=0 & vi[3] >=0){
    for(ki in c(k,30)){
      l <- ki+1
      for(i in 2:n){
        var[i]<-vi[1]+(vi[2]*gk(yc[i-1]/var[i-1],ki,l)+vi[3])*var[i-1]
      }
      ml <- mean(nfun(y2[2:n]-log(var[2:n]), div))
      if(ml<nml){
        nml <- ml
        if(ki==3){
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

freg <- function(x, a,b){

  u <- as.numeric(x>b)
  v <- as.numeric(x<a)

  c1 <- a-(2/(b-a)^3)*((-1/4)*a^4+(1/3)*(2*a+b)*a^3+(1/2)*(-a^2-2*a*b)*a^2+a^2*b*a)
  c2 <- (-1/(3*(b-a)^2))*(b-a)^3+b-(2/(b-a)^3)*((-1/4)*b^4+(1/3)*(2*a+b)*b^3+(1/2)*(-a^2-2*a*b)*b^2+a^2*b*b)

  g=x*v+(1-u-v)*((-1/(3*(b-a)^2))*(x-a)^3+x-(2/(b-a)^3)*((-1/4)*x^4+(1/3)*(2*a+b)*x^3+(1/2)*(-a^2-2*a*b)*x^2+a^2*b*x)+a-c1)+(c2+a-c1)*u

  g
}

gk <- function(x, k, l){

  bot <- 2*k-2*l
  a <- 1/bot
  b <- -2*l/bot
  c <- k^2/bot

  u <- as.numeric(x>l)
  v <- as.numeric(x<k)

  g<- x*v + (1-u-v)*(a*x^2 + b*x + c) + u*(a*l^2+b*l+c)

  g
}



