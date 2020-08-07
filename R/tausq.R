varin <- function(x)
{

  n <- length(x)
  a <- tausq(x)
  xsquare <- x^2
  int <- tausq(xsquare-1)
  i <- 1
  v <- x[1]^2

  while(v > a+int & i<30)
  {
    i <- i+1
    v <- x[i]^2
  }

  if(i==30)
  {
    error <- 1
  }

  v <- a
  v

}

tausq <- function(x)
{

  m <- length(x)

  s <- SEST2(x) # Global declaration.
  assign("Sestim", s, envir = .GlobalEnv)

  r <- RHO2(x/s)
  t <- mean(r)*s^2/0.4797

  t
}
