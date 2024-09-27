#' @title Robust Estimates for GARCH(1,1) Model
#'
#' @name robGarch
#' @aliases robGarch
#'
#' @description Methods for fitting a Garch(1,1) model with daily log return time series, using two methods of robust extended M-Estimates:
#' (1) maximizing modified likelihood function;
#' (2) M-Estimates with bounds to the propagation of the effect of one outlier on the subsequent predictors of the conditional variances.
#' returns a \code{robustGARCH} object
#'
#' @param data a time series of log returns, need to be numeric value.
#' @param methods robust M-Estimate method used for Garch(1,1) model, "M" and "BM", or non-robust M-Estimate method, "QML" and "MLE". Default is "BM".
#' @param fixed_pars a named numeric vector of parameters to be kept fixed during optimization, and they are needed for parameter estimation. For "M", the parameter should be c, which controls the modified loss function, user can use default c = .8; for "BM", the parameters should be c(c, k), where c is the same as in "M", user can use default c = 0.8,  and k (k > 0) is to control the robustness, the smaller k is, the more robust the method would be, user can use default k = 3.
#' @param optimizer optimizer used for optimization, one of "nloptr", "Rsolnp", "nlminb", default is "Rsolnp".
#' @param optimizer_x0 user-defined starting point for searching the optimum, c(x0_gamma, x0_alpha, x0_beta) or c(x0_gamma, x0_alpha, x0_beta, x0_shape) when the methods is "MLE". Default is "FALSE", where the starting point will be calculated instead of being user-defined.
#' @param optimizer_control list of control arguments passed to the optimizer. Default is list(trace=0). If wanting to print out the optimizer result, use list() instead.
#' @param stdErr_method method used to calculate standard error, one of "numDerive", "optim", "sandwich", default is "numDeriv" using hessian from numDeriv
#'
#' @return
#' A \code{robustGARCH} object(S3), the components of the object are:
#'     \item{data}{the log returns data object for the robustGARCH model to be fitted}
#'     \item{data_name}{the name of data variable input used}
#'     \item{methods}{the method called}
#'     \item{fixed_pars}{named numeric vector of fixed parameters used}
#'     \item{optimizer}{the optimizer called}
#'     \item{optimizer_x0}{user-defined or calculated starting point for searching the optimum}
#'     \item{optimizer_control}{the list of control arguments passed to the optimizer}
#'     \item{optimizer_result}{output of the called optimizer}
#'     \item{stdErr_method}{the method called to calculate standard error}
#'     \item{QML}{logical argument controlling the non-robustness of the fitting method}
#'     \item{fitted_pars}{Garch(1,1) parameter estimations output of the called method}
#'     \item{sigma}{the time series of the conditional standard deviation}
#'     \item{yt}{the time series of log(data^2)}
#'     \item{observed_I}{observed information matrix}
#'     \item{objective}{the optimal likihood value obtained by the optimizer}
#'     \item{time_elapsed}{the time used for the optimization routine}
#'     \item{message}{the message of the convergence status produced by the called solver}
#'     \item{standard_error}{standard erros of the fitted parameters using the method called}
#'     \item{t_value}{t-values of gamma, alpha, beta, shape as well for methods "MLE"}
#'     \item{p_value}{p-values of gamma, alpha, beta, shape as well for methods "MLE"}
#'
#' @details
#' The \code{robGarch} function fits a Garch(1, 1) model to a time series of log return data, using one of the two methods of robust extended M-Estimates with certain parameters specified by the user, with guidance and examples from the vignette. The user can also specify the optimizer used during optimization procesure, and the method used to calculate standard error for the fitted parameters.
#'
#' For details of the list of control arguments, please refer to \code{nloptr::nloptr}, \code{Rsolnp::solnp}, \code{nlminb}.
#'
#' @references Muler, Nora & Yohai, Victor. (2008). Robust estimates for GARCH models. Journal of Statistical Planning and Inference. 138. 2918-2940.
#'
#' @examples
#'
#'
#' data("gspc")
#' fit <- robGarch(gspc[1:604], methods="BM", fixed_pars = c(0.8, 3.0))
#'
#'
#' @rdname robustGARCH-robGarch
#' @export
# Garch(1,1) model fit function
robGarch <- function(data, methods = c("BM", "M", "QML", "MLE"), fixed_pars = c(0.8, 3.0), 
                     optimizer = c("Rsolnp", "nloptr", "nlminb"), optimizer_x0 = FALSE, 
                     optimizer_control = list(trace=0), stdErr_method = c("numDeriv", "optim", "sandwich")){
  # TODO: xts data takes a long time
  if(!is.numeric(data) || length(data)==0)
    stop("Data must be a numeric vector of non-zero length")

  methods = match.arg(methods)
  optimizer = match.arg(optimizer)
  stdErr_method = match.arg(stdErr_method)

  # Create a list to hold shared variables
  shared_vars <- list()

  shared_vars$methods <- methods

  if (methods == "QML" || methods == "MLE") {
    shared_vars$div <- 1.0
  } else {
    shared_vars$div <- fixed_pars[1]
  }

  if (methods == "BM") {
    shared_vars$k <- fixed_pars[2]
  } else {
    shared_vars$k <- 3.0
  }

  # Pass shared_vars to rgFit_local
  fit <- rgFit_local(data, optimizer, optimizer_x0, optimizer_control, shared_vars)

  # std_err calculation
  if(optimizer == "Rsolnp"){
    solution <- fit$optimizer_result$pars
    H <- fit$optimizer_result$hessian #[2:5, 2:5] for norm or [2:6, 2:6] for std.
  } else if(optimizer == "nloptr"){
    solution <- fit$optimizer_result$solution
    stop("use Rsolnp optimizer for now")
  } else if(optimizer == "nlminb"){
    solution <- fit$optimizer_result$par
    stop("use Rsolnp optimizer for now")
  }

  std_errors <- sqrt(diag(abs(solve(H)/length(data))))
  if(methods == "MLE"){
    standard_error <- c(std_errors[2:4], std_errors[6])
    fit$observed_I <- -H[c(2,3,4,6), c(2,3,4,6)]
  } else{
    standard_error <- std_errors[2:4]
    fit$observed_I <- -H[c(2,3,4), c(2,3,4)]
  }
  standard_error[1] <- standard_error[1] * (fit$fitted_pars[1] / solution[1])
  t_value <- fit$fitted_pars / standard_error
  p_value <- 2*(1-pnorm(abs(t_value)))
  
  ######## Change calculation method to above. ###########
  #if(stdErr_method == "numDeriv")
  #{
  #  stop("under development.")
  #  if(optimizer == "Rsolnp"){
  #    solution <- fit$optimizer_result$pars
  #  } else if(optimizer == "nloptr"){
  #    solution <- fit$optimizer_result$solution
  #  } else if(optimizer == "nlminb"){
  #    solution <- fit$optimizer_result$par
  #  }
  #  H <- numDeriv::hessian(Fnue, x = solution)/length(data)
  #  standard_error <- sqrt(diag(abs(solve(H)/length(data))))[1:3]
  #  t_value <- solution[1:3] / standard_error
  #  p_value <- 2*(1-pnorm(abs(t_value)))
  #}
  #if(stdErr_method == "optim")
  #{
  #  stop("under development")
  #  H <- optim(par = fit$optimizer_result$pars, fn = Fnue,
  #             method="L-BFGS-B",
  #             lower=c(-1, -1, -1, 0),
  #             upper=c(1, 1, 1, 2),
  #             hessian=TRUE)$hessian/length(data)
  #  standard_error <- sqrt(diag(abs(solve(H)/length(data))))[1:3]
  #  t_value <- fit$fitted_pars / standard_error
  #  p_value <- 2*(1-pnorm(abs(t_value)))
  #}
  #if(stdErr_method == "sandwich")
  #{
  #  stop("under development")
  #  Scores <- matrix(numDeriv::grad(Fnue, x = fit$optimizer_result$pars), nrow = 1)
  #  V <- (t(Scores) %*% Scores)/length(data)
  #  H <- optim(par = fit$optimizer_result$pars, fn = Fnue,
  #             method="L-BFGS-B",
  #             hessian=TRUE,
  #             control=list(maxit=1))$hessian/length(data)
  #  S <- solve(H) %*% V %*% solve(H)
  #  standard_error <- sqrt(diag(abs(S/length(data))))[1:3]
  #  t_value <- fit$fitted_pars / standard_error
  #  p_value <- 2*(1-pnorm(abs(t_value)))
  #}
  ##########################################

  if(methods == "MLE"){
    names(standard_error) <- c("gamma", "alpha", "beta", "shape")
    names(t_value) <- c("gamma", "alpha", "beta", "shape")
    names(p_value) <- c("gamma", "alpha", "beta", "shape")
  } else{
    names(standard_error) <- c("gamma", "alpha", "beta")
    names(t_value) <- c("gamma", "alpha", "beta")
    names(p_value) <- c("gamma", "alpha", "beta")
  }

  fit$data_name <- noquote(deparse(substitute(data)))
  fit$standard_error <- standard_error
  fit$t_value <- t_value
  fit$p_value <- p_value
  fit$fixed_pars <- fixed_pars

  structure(fit, class="robustGARCH")
}


robGarchDistribution <- function(param = c(8.76e-04, 0.135, 0.686), methods = c("BM", "M", "QML", "MLE"), fixed_pars = c(0.85, 3.0), optimizer = c("Rsolnp", "nloptr", "nlminb"), optimizer_x0 = FALSE, optimizer_control = list(), stdErr_method = c("numDeriv", "optim", "sandwich"), n = 2000, m = 100, rseed = 42){

  methods <- match.arg(methods)
  optimizer <- match.arg(optimizer)
  stdErr_method <- match.arg(stdErr_method)

  par(mfrow=c(2,2))
  spec <- rugarch::ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE))

  if(methods == "MLE"){
    if(length(param)!=4){stop("the parameters for std distribution should be gamma, alpha, beta, shape")}
    fixed <- param
    names(fixed) <- c("gamma", "alpha", "beta", "shape")
    fspec <- spec
    setfixed(fspec) <- fixed
    y <- rugarch::ugarchpath(fspec, n.sim = n, m.sim = m, rseed = 42)
    y. <- y@path$seriesSim

    qml_res <- matrix(0.0, nrow = m, ncol = 4)
    for( i in 1:m){
      y_ <- y.[((i-1)*n+1):(i*n)]
      fit <- robGarch(y_, methods = methods, fixed_pars = fixed_pars, optimizer=optimizer, optimizer_x0 = optimizer_x0, optimizer_control = optimizer_control, stdErr_method = stdErr_method)
      qml_res[i,1:4] <- fit$fitted_pars
    }

    d_omega <- density(qml_res[,1])
    d_omega
    plot(d_omega, main=paste("Parameter", expression( alpha ), "0 \n True value: ",fixed[1]), cex=.5)

    d_alpha1 <- density(qml_res[,2])
    d_alpha1
    plot(d_alpha1, main=paste("Parameter", expression( alpha ),"1 \n True value: ", fixed[2]), cex=.5)

    d_beta1 <- density(qml_res[,3])
    d_beta1
    plot(d_beta1, main=paste("Parameter", expression( beta ),"1 \n True value: ", fixed[3]), cex=.5)

    d_shape <- density(qml_res[,4])
    d_shape
    plot(d_shape, main=paste("Parameter shape\n True value: ", fixed[4]), cex=.5)
  } else{

    if(length(param)!=3){stop("the parameters for norm distribution should only be gamma, alpha, beta")}
    fixed <- param
    names(fixed) <- c("gamma", "alpha", "beta")
    fspec <- spec
    setfixed(fspec) <- fixed
    y <- rugarch::ugarchpath(fspec, n.sim = n, m.sim = m, rseed = 42)
    y. <- y@path$seriesSim

    qml_res <- matrix(0.0, nrow = m, ncol = 3)
    for( i in 1:m){
      y_ <- y.[((i-1)*n+1):(i*n)]
      fit <- robGarch(y_, methods = methods, fixed_pars = fixed_pars, optimizer=optimizer, optimizer_x0 = optimizer_x0, optimizer_control = optimizer_control, stdErr_method = stdErr_method)
      qml_res[i,1:3] <- fit$fitted_pars
    }

    d_omega <- density(qml_res[,1])
    d_omega
    plot(d_omega, main=paste("Parameter", expression( alpha ), "0 \n True value: ",fixed[1]), cex=.5)

    d_alpha1 <- density(qml_res[,2])
    d_alpha1
    plot(d_alpha1, main=paste("Parameter", expression( alpha ),"1 \n True value: ", fixed[2]), cex=.5)

    d_beta1 <- density(qml_res[,3])
    d_beta1
    plot(d_beta1, main=paste("Parameter", expression( beta ),"1 \n True value: ", fixed[3]), cex=.5)

  }
}
rgFit_local <- function(data, optimizer, optimizer_x0, optimizer_control, shared_vars){

  # unpack shared_vars
  methods <- shared_vars$methods

  start_time <- Sys.time()
  # optimizer/optimizer_control

  n <- length(data)
  v_data <- new_var(data, shared_vars)
  data_normalized <- data/sqrt(v_data)
  data_normalized <- data_normalized-median(data_normalized)
  for(i in 1:n){
    if(data_normalized[i] == 0){
      data_normalized[i] <- 10^(-10)
    }
  }
  vini <- new_var(data_normalized, shared_vars)

  res <- nEst(data_normalized, vini, optimizer, optimizer_x0, optimizer_control, shared_vars)
  optimizer_result <- res
    #if(optimizer == "fminsearch"){
      #stop("This feature is under active development, please use Rsolnp.")
      #fitted_pars <- res$fmin[1:3]
      #fitted_pars[1] <- fitted_pars[1]*v_data
      #names(fitted_pars) <- c("gamma", "alpha", "beta")
      #objective <- res$xopt
      #message <- NULL
    #} else
    if(optimizer == "Rsolnp"){
      if(methods == "MLE"){
        fitted_pars <- c(res$pars[1:3], res$pars[5])
        names(fitted_pars) <- c("gamma", "alpha", "beta", "shape")
      } else{
        fitted_pars <- res$pars[1:3]
        names(fitted_pars) <- c("gamma", "alpha", "beta")
      }
      fitted_pars[1] <- fitted_pars[1]*v_data
      objective <- res$values[length(res$values)]
      message <- res$convergence
    } else if (optimizer == "nloptr"){
      if(methods == "MLE"){
        fitted_pars <- c(res$solution[1:3], res$solution[5])
        names(fitted_pars) <- c("gamma", "alpha", "beta", "shape")
      } else{
        fitted_pars <- res$solution[1:3]
        names(fitted_pars) <- c("gamma", "alpha", "beta")
      }
      fitted_pars[1] <- fitted_pars[1]*v_data
      objective <- res$objective
      message <- res$message
    } else if (optimizer == "nlminb"){
      if(methods == "MLE"){
        fitted_pars <- c(res$par[1:3], res$par[5])
        names(fitted_pars) <- c("gamma", "alpha", "beta", "shape")
      } else{
        fitted_pars <- res$par[1:3]
        names(fitted_pars) <- c("gamma", "alpha", "beta")
      }
      fitted_pars[1] <- fitted_pars[1]*v_data
      objective <- res$objective
      message <- res$message
    } else{
      NA
    }
  sigma <- sigmaCal(fitted_pars, data, shared_vars)

  time_elapsed <- Sys.time() - start_time

  list(data=data,
       methods = methods,
       optimizer=optimizer,
       optimizer_x0=res$x0,
       optimizer_control=optimizer_control,
       optimizer_result=optimizer_result,
       fitted_pars = fitted_pars,
       objective=objective,
       time_elapsed=time_elapsed,
       message=message,
       sigma=sigma,
       yt=res$Muestram)
}
nEst <- function(y, vini, optimizer, optimizer_x0, optimizer_control, shared_vars){
  # unpack the shared_vars
  methods <- shared_vars$methods
  k <- shared_vars$k

  if(methods == "MLE"){
    std <- TRUE
  } else{
    std <- FALSE
  }

  n <- prod(length(y))
  yc <- y[1:n]^2
  y2 <- log(yc)

  shared_vars$Muestrac <- yc
  shared_vars$Muestram <- y2

  alfa0min <- 0.1
  alfa0max <- 1
  alfa1min <- 0.0
  alfa1max <- 1.
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
  if(std){

    if(length(optimizer_x0) > 1){
      x0 <- c(optimizer_x0[1:3], vini, optimizer_x0[4])
    } else{
      shapemin <- 3.0 # 1.0
      shapemax <- 31.
      nshape <- 30 # With 30 degree of freedom, the std would be very close to normal.
      # In the future, can seperate <30, >30 cases.
      lmshape <- (shapemax-shapemin)/nshape

      for(nj in 0:nalfa1){
        alfa1 <- alfa1min+nj*lmalfa1
        for(nk in 0:nbeta1){
          beta1 <- beta1min+nk*lmbeta1
          if(alfa1+beta1 <1 ){

            for(ni in 0:nalfa0){
              alfa0 <- alfa0min+ni*lmalfa0
              for(np in 0:nshape){
                shape <- shapemin+np*lmshape

                var <- rep(0.0, n)
                var[1] <- vini
                for(ki in c(k, 20)){
                  l <- ki+1
                  for(i in 2:n){
                    var[i] <- alfa0+(alfa1*rk(yc[i-1]/var[i-1],ki,l, shared_vars)+beta1)*var[i-1]
                  }

                  nml <- mean(nfun(y2[2:n]-log(var[2:n]), shared_vars, shape))
                  if(nml < ml){
                    flag <- 1
                    vi <- c(alfa0,alfa1,beta1,shape)
                    ml <- nml
                  }
                }
              }
            }
          }
        }
      }

      x0 <- c(vi[1:3], vini, vi[4])
    }

    lb <- c( 0., 0., 0., 0.0, 3.0) # 1.0 for lower bound for earlier version.
    ub <- c( 1.0, 1., 1., Inf, Inf)
  } else{

    if(length(optimizer_x0) > 1){
      x0 <- c(optimizer_x0[1:3], vini)
    } else{
      for(nj in 0:nalfa1){
        alfa1 <- alfa1min+nj*lmalfa1
        for(nk in 0:nbeta1){
          beta1 <- beta1min+nk*lmbeta1

          if(alfa1+beta1 <1 ){
            for(ni in 0:nalfa0){
              alfa0 <- alfa0min+ni*lmalfa0

              var <- rep(0.0, n)
              var[1] <- vini
              for(ki in c(k, 20)){
                l <- ki+1
                for(i in 2:n){
                  var[i] <- alfa0+(alfa1*rk(yc[i-1]/var[i-1],ki,l, shared_vars)+beta1)*var[i-1]
                }

                nml <- mean(nfun(y2[2:n]-log(var[2:n]), shared_vars))

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
      x0 <- c(vi, vini)
    }

    lb <- c( 0., 0., 0., 0.0)
    ub <- c( 1.0, 1.0, 1.0, Inf)
  }

  ### THIS IS FOR INVESGATING THE DIFF BETWEEN RUGARCH AND ROBUSTGARCH ###
  # print(x0)
  ### SHOULD BE DELETED ###

  if (optimizer == "nloptr"){

    #stop("This feature is under development, please use Rsolnp instead.")
    res <- nloptr::nloptr(x0 = x0,
                          eval_f = function(pars) Fnue(pars, shared_vars),
                          lb = lb,
                          ub = ub,
                          opts=optimizer_control)
    res$x0 <- x0
    res$Muestram <- y2

    return(res)

  } else if (optimizer == "Rsolnp"){

    res <- Rsolnp::solnp(pars = x0,
                         fun = function(pars) Fnue(pars, shared_vars),
                         LB = lb,
                         UB = ub,
                         ineqfun = function(vi){vi[2]+vi[3]},
                         ineqLB = 0.0,
                         ineqUB = 1.0,
                         control = optimizer_control)
    res$x0 <- x0
    res$Muestram <- y2

    return(res)

  } else if (optimizer == "nlminb"){

    res <- nlminb(start = x0,
                  objective = function(pars) Fnue(pars, shared_vars),
                  control = optimizer_control,
                  lower = lb,
                  upper = ub)
    res$x0 <- x0
    res$Muestram <- y2

    return(res)

  } else {
    NA
  }
}
new_var <- function(x, shared_vars){

  n <- length(x)
  tausq_x <- tau_sq(x, shared_vars)
  x_square <- x^2
  tausq_xsquare <- tau_sq(x_square-1, shared_vars)

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
tau_sq <- function(x, shared_vars){

  s <- s_est(x, shared_vars)
  r <- rho(x/s, shared_vars)
  t <- mean(r)*s^2/0.4797

  t
}
s_est <- function(x, shared_vars){

  b <- 1.625
  emed <- 0.675
  s <- 1 # Starting value
  eps <- 1.
  n <- 1

  m <- median(abs(x))/emed
  x <- x/m
  rho1 <- rho(x/0.405, shared_vars)
  a <- mean(rho1)/b
  v <- 1-a # Starting value
  si <- a # Starting value
  rho1 <- rho(x/(0.405*si), shared_vars) # Starting value
  a<- mean(rho1)/b # Starting value
  vi <- 1-a # Starting value
  AUX <- v * vi # Starting value

  while (eps> 0.005 & AUX > 0){
    n <- n+1
    s <- si
    v <- vi
    si <- a*s
    rho1 <- rho(x/(0.405 * si), shared_vars)
    a <- mean(rho1)/b
    vi <- 1-a
    AUX <- v*vi
    eps <- abs(s-si)/s
  }

  nsec <- 0
  while(eps>0.005){
    ns <- (s+si)/2
    rho1 <- rho(x/(0.405*ns), shared_vars)
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
rho <- function(x, shared_vars){

  # unpack shared_vars
  methods <- shared_vars$methods

  if(methods == "QML" || methods == "MLE" || methods == "M"){

    ps <- x^2/2

  } else{

    G1 <- -1.944
    G2 <- 1.728
    G3 <- -0.312
    G4 <- 0.016

    ax <- abs(x)
    u <- as.numeric(ax>3.0)
    v <- as.numeric(ax<2.0)
    w <- (1-u)*(1-v)
    ps <- v*x^2/2 + w*(G4*x^8/8 + G3*x^6/6 + G2*x^4/4 + G1*x^2/2 + 1.792) + 3.25*u
  }

  ps
}
nfun <- function(x, shared_vars, shape = 3.0){

  # unpack shared_vars
  methods <- shared_vars$methods

  # input here is w_t.
  b <- 4.3 #6.7428
  b1 <- 4.0 #b-0.5

  if(methods == "MLE"){
    # this is rho assuming z_t is std(shape).
    # this changes to MLE instead of QML.
    #x <- -log(gamma((shape+1)/2)/(sqrt((shape)*pi)*gamma(shape/2))) - x/2 + (shape+1)/2 *log(1+exp(x)/(shape))
    x <- -log(gamma((shape+1)/2)/(sqrt(shape-2)*gamma(shape/2))) - x/2 + (shape+1)/2 *log(1+exp(x)/(shape-2))
  } else{
    # assume z_t is standard norm.
    x <- (exp(x)-x +log(2*pi))/2
  }

  ps <- freg(x, b1, b, shared_vars)

  ps
}
sigmaCal <- function(pars, data, shared_vars){

  # unpack the shared_vars
  k <- shared_vars$k

  n <- prod(length(data))

  var <- rep(0.0, n)
  var[1] <- pars[1]/(1-pars[3])

  if(pars[1] >0 & pars[2] >=0 & pars[3] >=0){
    l <- k+1
    for(i in 2:n){
      # h_c(t) in the paper.
      var[i]<-pars[1]+(pars[2]*rk(data[i-1]^2/var[i-1],k,l, shared_vars)+pars[3])*var[i-1]
    }
  }

  var
}
Fnue <- function(start_pars, shared_vars){

  # unpack the shared_vars
  methods <- shared_vars$methods
  k <- shared_vars$k

  y2 <- shared_vars$Muestram
  yc <- shared_vars$Muestrac

  n <- prod(length(y2))

  vi <- start_pars[1:3]
  vini <- start_pars[4]
  if(methods == "MLE"){
    shape <- start_pars[5]
  }

  var <- rep(0.0, n)
  var[1] <- vini

  nml <- 10^7

  if(vi[1] >0 & vi[2] >=0 & vi[3] >=0){
    for(ki in c(k,20)){

      l <- ki+1
      for(i in 2:n){
        var[i]<-vi[1]+(vi[2]*rk(yc[i-1]/var[i-1],ki,l, shared_vars)+vi[3])*var[i-1]
      }
      if(methods == "MLE"){
        ml <- mean(nfun(y2[2:n]-log(var[2:n]), shared_vars, shape))
      } else{
        ml <- mean(nfun(y2[2:n]-log(var[2:n]), shared_vars))
      }

      if(is.nan(ml) || is.nan(nml)){
        break
      }
      else if(ml<nml){ nml <- ml}
    }
  }

  nml
}
freg <- function(x, a, b, shared_vars){
  # the rho function.

  # unpack shared_vars
  methods <- shared_vars$methods
  div <- shared_vars$div

  if(methods == "QML" || methods == "MLE" || a == b){
    # used to be return exp(w/div)-w/div, now correct it to be the following as stated in the paper.
    g <- x

  } else{

    x <- x/div

    u <- as.numeric(x>b)
    v <- as.numeric(x<a)

    ba <- b-a
    c1 <- a-(2*a^3/ba^3)*(b/3-a/12) #a-(2/(b-a)^3)*((-1/4)*a^4+(1/3)*(2*a+b)*a^3+(1/2)*(-a^2-2*a*b)*a^2+a^2*b*a)
    c2 <- -ba/3 +b -(2*b^2/ba^3)*(b^2/12-a*b/3+a^2/2) #(-1/(3*(b-a)^2))*(b-a)^3+b-(2/(b-a)^3)*((-1/4)*b^4+(1/3)*(2*a+b)*b^3+(1/2)*(-a^2-2*a*b)*b^2+a^2*b*b)

    g <- v*x +
      (1-u-v)*(x -(x-a)^3/(3*ba^2) -(2*a^2*b*(x-a))/ba^3 +2*((x^4-a^4)/4 -(2*a+b)*(x^3-a^3)/3 +(a^2+2*a*b)*(x^2-a^2)/2)/ba^3) +
      u*(c2+a-c1)

    g <- g*div
  }

  g
}
rk <- function(x, ki, l, shared_vars){

  # unpack shared_vars
  methods <- shared_vars$methods
  # k should not be used as a global variable as it has conflict with ki (previously k) in the function.
  # k <- shared_vars$k

  if(methods == "QML" || methods == "MLE" || methods == "M" || ki == l){

    g <- x

  } else{
    bot <- 2*ki-2*l
    a <- 1/bot
    b <- -2*l/bot
    c <- ki^2/bot

    u <- as.numeric(x>l)
    v <- as.numeric(x<ki)

    g<- x*v + (1-u-v)*(a*x^2 + b*x + c) + u*(a*l^2+b*l+c)
  }

  g
}



