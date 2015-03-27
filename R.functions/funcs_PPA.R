# 3 different intensity functions to model the start and end 
# times as a non-homogenuous Poisson Process
lambda1 <- function(t, theta){
  theta[1] + 0*t
}
lambda2 <- function(t, theta) {
  theta[1] + theta[2] * sin(2*pi*t) + theta[3] * cos(2*pi*t)
} 
lambda3 <- function(t, theta) { 
  theta[1] + theta[2] * dnorm(tan(pi*t+theta[3]), mean = theta[4], sd = theta[5])
}

# fix to integrate taking into acount the jump between years in ss.frac and se.frac
library(stats4)
my.integrate <- function(lambda, theta, lower, upper){
  if(upper-lower>0){
    I <- integrate(lambda, theta, lower=lower,upper=upper)[[1]]
  }else{
    I <- integrate(lambda, theta, lower = upper, upper = 1)[[1]] + 
         integrate(lambda, theta, lower = 0, upper = lower)[[1]]
  }
  I
}

# maximum likelihood fit of gaps to a pp.distribution with intensity 
# function lambda1
ll.lam1 <- function(theta){
  print(theta)
  l   <- length(se.frac)
  vI <- vector("numeric")
  vlg <- vector("numeric")
  for (i in 1:(l-1)){
    vI[i]  <- my.integrate(lambda1,theta,lower=se.frac[i],upper=ss.frac[i+1])
    vlg[i] <- log(lambda1(ss.frac[i+1],theta))
  }
  loglike <-  sum(vI) - sum(vlg)
  loglike
} 

# maximum likelihood fit of gaps to a pp.distribution with intensity 
# function lambda2
llike.lam2 <- function(theta1, theta2, theta3){
  theta <- c(theta1, theta2, theta3)
  print(theta)
  l   <- length(se.frac)
  vI <- vector("numeric")
  vlg <- vector("numeric")
  for (i in 1:(l-1)){
    vI[i]  <- my.integrate(lambda2,theta,lower=se.frac[i],upper=ss.frac[i+1])
    vlg[i] <- log(lambda2(ss.frac[i+1],theta))
  }
  loglike <-  sum(vI) - sum(vlg)
  loglike
}

# maximum likelihood fit of gaps to a pp.distribution with intensity 
# function lambda3. In order to converge is necessary to estimate 3
# parameters first, and then add the other 2. 
lam3.aux <- function(t, theta) {
  theta[1] + theta[2] * dnorm(tan(pi*t+theta[3]), mean = 0, sd = 1)
}

llike.lam3.aux <- function(theta1, theta2, theta3){
  theta <- c(theta1, theta2, theta3)
  print(theta)
  l   <- length(se.frac)
  vI <- vector("numeric")
  vlg <- vector("numeric")
  for (i in 1:(l-1)){
    vI[i]  <- my.integrate(lam3.aux,theta,lower=se.frac[i],upper=ss.frac[i+1])
    vlg[i] <- log(lam3.aux(ss.frac[i+1],theta))
  }
  loglike <-  sum(vI) - sum(vlg)
  loglike
}
llike.lam3 <- function(theta1, theta2, theta3, theta4, theta5){
  theta <- c(theta1, theta2, theta3, theta4, theta5)
  print(theta)
  l   <- length(se.frac)
  vI <- vector("numeric")
  vlg <- vector("numeric")
  for (i in 1:(l-1)){
    vI[i]  <- my.integrate(lambda3,theta,lower=se.frac[i],upper=ss.frac[i+1])
    vlg[i] <- log(lambda3(ss.frac[i+1],theta))
  }
  loglike <-  sum(vI) - sum(vlg)
  loglike
}

p.pp.cond <- function(x, lambda, theta, ts){
# Probability function of a non-homogenous Possion Process with density lambda
# (computes F(ts+x|ts) where
# F(ts+x|ts) = 1 - \exp(- \int_ts^{ts+x} \lambda(s |theta) ds) )
  if(x<0){ stop("negative integration time")
  }else if(x<1e-14){ I <- 0
  }else{  I <- my.integrate(lambda,theta,lower=ts,upper=ts+x)}
  p  <- 1-exp(-I)
  p
}

d.pp.cond <- function(x, lambda, theta, ts){
# Density function of a non-homogenous Possion Process with density lambda
# (computes F'(ts+x|ts) where
# F(ts+x|ts) = 1 - \exp(- \int_ts^{ts+x} \lambda(s |theta) ds) )
  if(x<0){ stop("negative integration time")
  }else if(x<1e-14){ I <- 0
  }else{  I <- my.integrate(lambda,theta,lower=ts,upper=ts+x)}
  d  <- exp(-I)*lambda(ts+x,theta)
}

bisection <- function( f, a, b, tol=1e-5){
# bisection method to find a zero of the function  fx between a, b

  fa    <- f(a)
  fb    <- f(b)
  while(abs(f(a))> tol) {
    fhalf <- f((a+b)/2)
    if(fa*fb>0){
      warning(c("bisecion: a=", a," b=", b, " f(a)=", f(a)," f(b)=",f(b)))
      return(NaN)
    }
    if(fa*fhalf <0){
      b <- (a+b)/2; fb <- fhalf
    }else{
      a <- (a+b)/2; fa <- fhalf
    }
  }
  (a+b)/2
}

r.pp.cond <- function(lambda, theta, ts, tol=1e-12){
# Random number generator distributed as a non-homogenous Possion Process 
# with density lambda (computes x ~ F(ts+x|ts) where
# F(ts+x|ts) = 1 - \exp(- \int_ts^{ts+x} \lambda(s |theta) ds) )
  u1 = runif(1)
  H  <- function(x) (p.pp.cond(x, lambda, theta, ts) - u1)
  Hp <- function(x) d.pp.cond(x, lambda, theta, ts)
  g0 <- bisection(H, 0, 2, tol)
  g0
}


test.r.pp.cond <- function(lambda, thet.lam, ts, nsamp = 1000){
#  s.pp.cond <- vector("numeric")
  s.pp.cond <- vector("numeric")
  for (i in 1:nsamp){
    s.pp.cond[i] <- r.pp.cond(lambda,thet.lam,ts)
    print(c(i,s.pp.cond[i]))
  }
  hist(s.pp.cond, prob=T)
  xx <- seq(0,1,0.0001)
  yy <- sapply(xx, function(x) {d.pp.cond(x,lambda,thet.lam,ts)})
  lines(xx,yy)

  s.pp.cond
}

#dat.test <- test(lambda1,thet.lam1,0,nsamp=5000)
#dat.test <- test(lambda2,thet.lam2,0,nsamp=1000)
#dat.test <- test(lambda3,thet.lam3,0,nsamp=1000)
#str(dat.test)
#which(is.na(dat.test))
#dat.test[which(is.na(dat.test)),]

