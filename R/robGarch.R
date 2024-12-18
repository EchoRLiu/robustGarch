#' @importFrom stats median pnorm density nlminb
#' @importFrom graphics par
#' @title Robust GARCH(1,1) Model Estimation
#'
#' @name robGarch
#'
#' @description Computes "BM" robust Garch(1,1) model parameter
#' estimate by using a bounded objective function and a bounded
#' conditional variance recursion.  Alternatively, it computes:
#' (1) "M" estimates by using only the bounded objective function,
#' (2) "QML" estimates based on a typically incorrect assumption
#' of normally distributed innovations, (3) "t-MLE" estimates based
#' on an assumption of an innovations t-distributed MLE with unknown
#' location, scale,and degrees of freedom parameters.
#' CHECK IF (3) IS CORRECT.
#'
#' @references Muler, N. and Yohai, V. (2008). Robust estimates
#' for GARCH models. Journal of Statistical Planning and Inference,
#' 138, 2918-2940.
#'
#' @param data an xts object
#' @param fitMethod character valued name of fitting method,
#' one of "BM", "M" "QML" or "tMLE", with "BM" the default value.
#' @param robTunePars a numeric vector c(cM,cFlt) that controls
#' the extent of fitMethod robustness, with default c(0.8,3.0).
#' @param optChoice character valued optChoice name, one of
#' "Rsolnp", "nloptr", "nlminb", with default "Rsolnp".
#' @param initialPars numeric user-defined initial parameters
#' c(gamma0, alpha0, beta0) for use by optChoice, with default
#' values c(0.0005, 0.15, 0.75).
#' @param optControl list of arguments passed to optChoice, with
#' default \code{list(trace=0)}.
#' @param SEmethod character valued name of standard error method,
#' one of "numDeriv", "optim", "sandwich", with default "numDeriv".
#'
#' @details The "BM" fit method delivers the highest robustness by
#' using a half-Huber psi function to bound the normal distribution
#' log-likelihood, and using a Huber psi function to prevent the
#' propagation of influential outliers in the variance recursion.
#' The "M" method is obtained by dropping the BM bounding of the
#' variance recursion, and is therefore less robust toward outliers.
#'
#' ECHO OR DAN, PLEASE PROVIDE DETAILS FOR optControl.
#' For details of the list of control arguments, please refer to
#' \code{nloptr::nloptr}, \code{Rsolnp::solnp}, \code{nlminb}.
#' The SEmethod default "numDeriv" is based on the Hessian from the
#' optimization.
#'
#' @return
#' A list object of class \dQuote{robustGarch} with components:
#' \item{data}{the input xts object}
#' \item{fitMethod}{the the fitMethod specified}
#' \item{robtunePars}{the robtunePars specified}
#' \item{initialPars}{the initialPars specified}
#' \item{optChoice}{the optChoice specified}
#' \item{coefEstimates}{computed parameter estimates}
#' \item{sigma}{conditional standard deviation xts class time series}
#' \item{SEmethod}{the specidied of calculating standard errors}
#' \item{observedInfoMat}{observed information matrix}
#' \item{optDetails}{a list containing the optChoice specified,
#' the control values specified, and the optChoice minimized
#' objective, and convergence status message}
#'
#' @rdname robustGARCH-robGarch
#' @export
#'
#' @examples
#' data("gspc")
#' fit <- robGarch(gspc[1:604], fitMethod = "BM")
#' summary(fit)
#'
robGarch <- function(data,
                     fitMethod = c("BM", "M", "QML", "MLE"),
                     robTunePars = c(0.8, 3.0),
                     optChoice = c("Rsolnp", "nloptr", "nlminb"),
                     initialPars = c(0.0005, 0.15, 0.75),
                     SEmethod = c("numDeriv", "optim", "sandwich"),
                     optControl = list(trace=0))
  {

  # if(!is.numeric(data) || length(data)==0)
  #   stop("Data must be a numeric vector of non-zero length")

  # the calculation will use data_, which only has the return data and without dates
  # but we will keep data to retain the original dates
  data_ = zoo::coredata(data)

  fitMethod = match.arg(fitMethod)
  optChoice = match.arg(optChoice)
  SEmethod = match.arg(SEmethod)
  print(SEmethod)

  # Create a list to hold shared variables
  shared_vars <- list()

  shared_vars$fitMethod <- fitMethod

  if (fitMethod == "QML" || fitMethod == "MLE") {
    shared_vars$div <- 1.0
  } else {
    shared_vars$div <- robTunePars[1]
  }

  if (fitMethod == "BM") {
    shared_vars$k <- robTunePars[2]
  } else {
    shared_vars$k <- 3.0
  }

  # Pass shared_vars to rgFit_local
  fit <- rgFit_local(data, optChoice, initialPars, optControl, shared_vars)

  # std_err calculation
  if(optChoice == "Rsolnp"){
    solution <- fit$optChoice_result$pars
    H <- fit$optChoice_result$hessian #[2:5, 2:5] for norm or [2:6, 2:6] for std.
  } else if(optChoice == "nloptr"){
    solution <- fit$optChoice_result$solution
    stop("use Rsolnp optChoice for now")
  } else if(optChoice == "nlminb"){
    solution <- fit$optChoice_result$par
    stop("use Rsolnp optChoice for now")
  }

  std_errors <- sqrt(diag(abs(solve(H)/length(data_))))
  if(fitMethod == "MLE"){
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
  #if(SEmethod == "numDeriv")
  #{
  #  stop("under development.")
  #  if(optChoice == "Rsolnp"){
  #    solution <- fit$optChoice_result$pars
  #  } else if(optChoice == "nloptr"){
  #    solution <- fit$optChoice_result$solution
  #  } else if(optChoice == "nlminb"){
  #    solution <- fit$optChoice_result$par
  #  }
  #  H <- numDeriv::hessian(Fnue, x = solution)/length(data)
  #  standard_error <- sqrt(diag(abs(solve(H)/length(data))))[1:3]
  #  t_value <- solution[1:3] / standard_error
  #  p_value <- 2*(1-pnorm(abs(t_value)))
  #}
  #if(SEmethod == "optim")
  #{
  #  stop("under development")
  #  H <- optim(par = fit$optChoice_result$pars, fn = Fnue,
  #             method="L-BFGS-B",
  #             lower=c(-1, -1, -1, 0),
  #             upper=c(1, 1, 1, 2),
  #             hessian=TRUE)$hessian/length(data)
  #  standard_error <- sqrt(diag(abs(solve(H)/length(data))))[1:3]
  #  t_value <- fit$fitted_pars / standard_error
  #  p_value <- 2*(1-pnorm(abs(t_value)))
  #}
  #if(SEmethod == "sandwich")
  #{
  #  stop("under development")
  #  Scores <- matrix(numDeriv::grad(Fnue, x = fit$optChoice_result$pars), nrow = 1)
  #  V <- (t(Scores) %*% Scores)/length(data)
  #  H <- optim(par = fit$optChoice_result$pars, fn = Fnue,
  #             method="L-BFGS-B",
  #             hessian=TRUE,
  #             control=list(maxit=1))$hessian/length(data)
  #  S <- solve(H) %*% V %*% solve(H)
  #  standard_error <- sqrt(diag(abs(S/length(data))))[1:3]
  #  t_value <- fit$fitted_pars / standard_error
  #  p_value <- 2*(1-pnorm(abs(t_value)))
  #}
  ##########################################

  if(fitMethod == "MLE"){
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
  fit$robTunePars <- robTunePars
  fit$fitMethod <- fitMethod

  structure(fit, class="robustGARCH")
}


robGarchDistribution <- function(
  param = c(8.76e-04, 0.135, 0.686),
  fitMethod = c("BM", "M", "QML", "MLE"),
  robTunePars = c(0.85, 3.0),
  optChoice = c("Rsolnp", "nloptr", "nlminb"),
  initialPars = FALSE,
  optControl = list(),
  SEmethod = c("numDeriv", "optim", "sandwich"),
  n = 2000,
  m = 100,
  rseed = 42) {

  fitMethod <- match.arg(fitMethod)
  optChoice <- match.arg(optChoice)
  SEmethod <- match.arg(SEmethod)

  par(mfrow=c(2,2))
  spec <- rugarch::ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE))

  if(fitMethod == "MLE"){
    if(length(param)!=4){stop("the parameters for std distribution should be gamma, alpha, beta, shape")}
    fixed <- param
    names(fixed) <- c("gamma", "alpha", "beta", "shape")
    fspec <- spec
    rugarch::setfixed(fspec) <- fixed
    y <- rugarch::ugarchpath(fspec, n.sim = n, m.sim = m, rseed = 42)
    y. <- y@path$seriesSim

    qml_res <- matrix(0.0, nrow = m, ncol = 4)
    for( i in 1:m){
      y_ <- y.[((i-1)*n+1):(i*n)]
      fit <- robGarch(
        y_,
        fitMethod = fitMethod,
        robTunePars = robTunePars,
        optChoice = optChoice,
        initialPars = initialPars,
        optControl = optControl,
        SEmethod = SEmethod
      )
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
    rugarch::setfixed(fspec) <- fixed
    y <- rugarch::ugarchpath(fspec, n.sim = n, m.sim = m, rseed = 42)
    y. <- y@path$seriesSim

    qml_res <- matrix(0.0, nrow = m, ncol = 3)
    for( i in 1:m){
      y_ <- y.[((i-1)*n+1):(i*n)]
      fit <- robGarch(y_, fitMethod = fitMethod, robTunePars = robTunePars, optChoice=optChoice, initialPars = initialPars, optControl = optControl, SEmethod = SEmethod)
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
rgFit_local <- function(data, optChoice, initialPars, optControl, shared_vars){

  # unpack shared_vars
  fitMethod <- shared_vars$fitMethod

  # for minimum code change, data_ and data naming are exchanged here
  data_ = data # retain the dates and structure
  data = zoo::coredata(data) # for calculation

  start_time <- Sys.time()
  # optChoice/optControl

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

  res <- nEst(data_normalized, vini, optChoice, initialPars, optControl, shared_vars)
  optChoice_result <- res
    #if(optChoice == "fminsearch"){
      #stop("This feature is under active development, please use Rsolnp.")
      #fitted_pars <- res$fmin[1:3]
      #fitted_pars[1] <- fitted_pars[1]*v_data
      #names(fitted_pars) <- c("gamma", "alpha", "beta")
      #objective <- res$xopt
      #message <- NULL
    #} else
    if(optChoice == "Rsolnp"){
      if(fitMethod == "MLE"){
        fitted_pars <- c(res$pars[1:3], res$pars[5])
        names(fitted_pars) <- c("gamma", "alpha", "beta", "shape")
      } else{
        fitted_pars <- res$pars[1:3]
        names(fitted_pars) <- c("gamma", "alpha", "beta")
      }
      fitted_pars[1] <- fitted_pars[1]*v_data
      objective <- res$values[length(res$values)]
      message <- res$convergence
    } else if (optChoice == "nloptr"){
      if(fitMethod == "MLE"){
        fitted_pars <- c(res$solution[1:3], res$solution[5])
        names(fitted_pars) <- c("gamma", "alpha", "beta", "shape")
      } else{
        fitted_pars <- res$solution[1:3]
        names(fitted_pars) <- c("gamma", "alpha", "beta")
      }
      fitted_pars[1] <- fitted_pars[1]*v_data
      objective <- res$objective
      message <- res$message
    } else if (optChoice == "nlminb"){
      if(fitMethod == "MLE"){
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

  list(data=data_,
       fitMethod = fitMethod,
       optChoice=optChoice,
       initialPars=res$x0,
       optControl=optControl,
       optChoice_result=optChoice_result,
       fitted_pars = fitted_pars,
       objective=objective,
       time_elapsed=time_elapsed,
       message=message,
       sigma=sigma,
       yt=res$Muestram)
}
nEst <- function(y, vini, optChoice, initialPars, optControl, shared_vars){
  # unpack the shared_vars
  fitMethod <- shared_vars$fitMethod
  k <- shared_vars$k

  if(fitMethod == "MLE"){
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

    if(length(initialPars) > 1){
      x0 <- c(initialPars[1:3], vini, initialPars[4])
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

    if(length(initialPars) > 1){
      x0 <- c(initialPars[1:3], vini)
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

  if (optChoice == "nloptr"){

    #stop("This feature is under development, please use Rsolnp instead.")
    res <- nloptr::nloptr(x0 = x0,
                          eval_f = function(pars) Fnue(pars, shared_vars),
                          lb = lb,
                          ub = ub,
                          opts=optControl)
    res$x0 <- x0
    res$Muestram <- y2

    return(res)

  } else if (optChoice == "Rsolnp"){

    res <- Rsolnp::solnp(pars = x0,
                         fun = function(pars) Fnue(pars, shared_vars),
                         LB = lb,
                         UB = ub,
                         ineqfun = function(vi){vi[2]+vi[3]},
                         ineqLB = 0.0,
                         ineqUB = 1.0,
                         control = optControl)
    res$x0 <- x0
    res$Muestram <- y2

    return(res)

  } else if (optChoice == "nlminb"){

    res <- nlminb(start = x0,
                  objective = function(pars) Fnue(pars, shared_vars),
                  control = optControl,
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
  fitMethod <- shared_vars$fitMethod

  if(fitMethod == "QML" || fitMethod == "MLE" || fitMethod == "M"){

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
  fitMethod <- shared_vars$fitMethod

  # input here is w_t.
  b <- 4.3 #6.7428
  b1 <- 4.0 #b-0.5

  if(fitMethod == "MLE"){
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
  fitMethod <- shared_vars$fitMethod
  k <- shared_vars$k

  y2 <- shared_vars$Muestram
  yc <- shared_vars$Muestrac

  n <- prod(length(y2))

  vi <- start_pars[1:3]
  vini <- start_pars[4]
  if(fitMethod == "MLE"){
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
      if(fitMethod == "MLE"){
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
  fitMethod <- shared_vars$fitMethod
  div <- shared_vars$div

  if(fitMethod == "QML" || fitMethod == "MLE" || a == b){
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
  fitMethod <- shared_vars$fitMethod
  # k should not be used as a global variable as it has conflict with ki (previously k) in the function.
  # k <- shared_vars$k

  if(fitMethod == "QML" || fitMethod == "MLE" || fitMethod == "M" || ki == l){

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



