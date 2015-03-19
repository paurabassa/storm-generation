# Read data
waves <- read.csv("./Clean.Data/waves-3hours.csv")
waves$Date <- as.POSIXct(strptime(waves$Date, "%Y-%m-%d %H:%M:%S"))

# Extract the point process of excedances (storms)
source("functions_pots.R")
pots.storm <- ExtractPOTPointProcess(waves$hs,2.5)
storm.start <- as.POSIXct(waves$Date[pots.storm$p.exc])
storm.end   <- as.POSIXct(waves$Date[pots.storm$p.exc + pots.storm$c.siz])

# Auxiliary function to compute the storm start and ends as fraction of the year
yearfraction <- function(date){
   print(date) 
   aux   <- as.POSIXlt(date)
   year  <- aux$year
   lyear <- difftime(as.POSIXct(strptime(paste(1901+year,"-01-01",sep=""),"%Y-%m-%d")),  
                     as.POSIXct(strptime(paste(1900+year,"-01-01",sep=""),"%Y-%m-%d")),
                     units="hours")
   fyear <- (24*aux$yday + aux$hour)/lyear[[1]]
}

# Computation of the year fraction times and sanity test
ss.frac <- sapply(storm.start,yearfraction)
se.frac <- sapply(storm.end  ,yearfraction)
l <- length(ss.frac)
which(ss.frac[2:l]-se.frac[1:(l-1)]<0)

# histogram of diferences taking into acount the jump from one year to the next
faux<- function(x){ if(x<0){ y <-x+1}else{ y<-x}; y}
gaps <- sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux)
hist(gaps,breaks=50)
nyears <- as.POSIXlt(tail(waves$Date,n=1))$year - as.POSIXlt(waves$Date[1])$year

# check if there is a correlation between storm duration (c.siz) and the 
# inter arrival times
cor(sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux),pots.storm$c.siz[1:(l-1)])

# intensity functions to model the start and end times as a 
# non-homogenuous Poisson Process
lambda1 <- function(t, theta){
  theta[1] + 0*t
}
lambda2 <- function(t, theta) {
  theta[1] + theta[2] * sin(2*pi*t) + theta[3] * cos(2*pi*t)
} 
lambda3 <- function(t, theta) { 
  theta[1] + theta[2] * dnorm(tan(pi*t+theta[3]), mean = theta[4], sd = theta[5])
}

# plot of the intensity functions above
x <- seq(-1,1,0.01)
y1 <- sapply(x,function(x) { lambda1(x,1)}) 
y2 <- sapply(x,function(x) { lambda2(x,c(1,0.5,0))})
y3 <- sapply(x,function(x) { lambda3(x,c(0,1,0,0,1))})
plot(x,y1,type="l",clim=c(0,2))
lines(x,y2,lty=2)
lines(x,y3,lty=3)

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
fit.lam1  <- mle(ll.lam1,list(theta=1))
thet.lam1 <- coef(fit.lam1)[[1]]


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
fit.lam2  <- mle( llike.lam2, start = list(theta1 = 1, theta2=0.4, theta3=0.4))
thet.lam2 <- coef(fit.lam2)


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
llike.lam3.aux(0,1,0)
fit.lam3.aux <- mle( llike.lam3.aux, method="L-BFGS-B", lower=c(0.01,0.01,0),
                  start = list(theta1 = 0, theta2=1, theta3=0.))
thet.lam3.aux <- coef(fit.lam3.aux)

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
llike.lam3(0,1,0,0,1)
fit.lam3  <- mle( llike.lam3, method="BFGS", lower=c(0,0.01,0,-Inf,0),
                  start = list(theta1 = thet.lam3.aux[[1]], theta2=thet.lam3.aux[[2]], 
                               theta3 = thet.lam3.aux[[3]], theta4 = 0, theta5=1))
thet.lam3 <- coef(fit.lam3)

#Comparasion of the 3 fits
x      <- seq(0,1,0.001)
y1     <- sapply(x,function(x) { lambda1(x,thet.lam1)})
area   <- sum(y1)/length(y1)
y1     <- y1/area

y2     <- sapply(x,function(x) { lambda2(x,thet.lam2)})
area   <- sum(y2)/length(y2)
y2     <- y2/area

y3     <- sapply(x,function(x) { lambda3(x,thet.lam3)})
area   <- sum(y3)/length(y3)
y3     <- y3/area

hist(ss.frac, prob = T)
lines(x, y1, lty=1)
lines(x, y2, lty=2)
lines(x, y3, lty=3)


# preparation of the quantile-quantile plots

# auxiliary functions
xIpd.pp.cond <- function(x, lambda, theta, ts){ 
# Density and propability function of a non-homogenous Possion
# Function to compute F(ts+x|ts) and F'(ts+x|ts) where
# F(ts+x|ts) = 1 - \exp(- \int_ts^{ts+x} \lambda(s |theta) ds
  if(x<0){ stop("negative integration time")
  }else if(x<1e-14){ I <- 0
  }else{  I <- my.integrate(lambda,theta,lower=ts,upper=ts+x)}
  p  <- 1-exp(-I)
  d  <- exp(-I)*lambda(ts+x,theta)
  c(x,I,p,d)
}

p.pp.cond <- function(x, lambda, theta, ts){
  if(x<0){ stop("negative integration time")
  }else if(x<1e-14){ I <- 0
  }else{  I <- my.integrate(lambda,theta,lower=ts,upper=ts+x)}
  p  <- 1-exp(-I)
  p
}

d.pp.cond <- function(x, lambda, theta, ts){
  if(x<0){ stop("negative integration time")
  }else if(x<1e-14){ I <- 0
  }else{  I <- my.integrate(lambda,theta,lower=ts,upper=ts+x)}
  d  <- exp(-I)*lambda(ts+x,theta)
  d
}

xIpd.pp.cond(0.0001,lambda1, thet.lam1, 0.000)
m <- sapply(seq(0.001,1,0.001),function(x){ xIpd.pp.cond(x,lambda3, thet.lam3, 0.400)})
plot(m[1,],m[2,])
plot(m[1,],log(1-m[3,]))
plot(m[1,],log(m[4,]))

m <- sapply(seq(0.001,1,0.001),function(x){ xIpd.pp.cond(x,lambda2, thet.lam2, 0.000)})
plot(m[1,],m[3,])
m <- sapply(seq(0.001,1,0.001),function(x){ xIpd.pp.cond(x,lambda2, thet.lam2, 0.5)})
plot(m[1,],m[3,])

p.pp.marginal <- function(x, lambda, theta){ 
# probability integrated over all possible ts \in [0,1] 
  h<- 0.0001
  nodes<- sapply(seq(0,1,h), function(ts){ p.pp.cond(x,lambda,theta,ts)})
  sum(nodes)*h
}

d.pp.marginal <- function(x, lambda, theta){ 
# density integrated over all possible ts \in [0,1] 
  h<- 0.0001
  nodes<- sapply(seq(0,1,h), function(ts){ d.pp.cond(x,lambda,theta,ts)})
  (sum(nodes)-nodes[1])*h 
}

prob <- sapply(seq(0,1,0.01),function(t){p.pp.marginal(t,lambda1,thet.lam1)}) 
dens  <- sapply(seq(0,1,0.01),function(t){d.pp.marginal(t,lambda1,thet.lam1)})
prob2 <- sapply(seq(0,1,0.01),function(t){p.pp.marginal(t,lambda2,thet.lam2)})
dens2 <- sapply(seq(0,1,0.01),function(t){d.pp.marginal(t,lambda2,thet.lam2)})
prob3 <- sapply(seq(0,1,0.01),function(t){p.pp.marginal(t,lambda3,thet.lam3)})
dens3 <- sapply(seq(0,1,0.01),function(t){d.pp.marginal(t,lambda3,thet.lam3)})

hist(gaps,prob=T)
lines(seq(0,1,0.01),dens,  col="green")
lines(seq(0,1,0.01),dens2, col="blue")
lines(seq(0,1,0.01),dens3, col="red")

ecdf <- ecdf(gaps)
plot(ecdf)
lines(seq(0,1,0.01),prob,  col="green")
lines(seq(0,1,0.01),prob2, col="blue")
lines(seq(0,1,0.01),prob3, col="red")

