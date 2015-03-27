# 1st auxiliary functions 
Trapez <- function( x, f){
  area <- 0
  for (i in 1:(length(x)-1))
    area<- area + (x[i+1] - x[i])*(f[i+1]+f[i])/2
  area
}

CDF <- function(x,f){
  F<- numeric(length(x))
  F[1] <- 0
  for (i in 1:(length(x)-1))
    F[i+1] <- F[i] + (x[i+1] - x[i])*(f[i+1]+f[i])/2
  F
}


