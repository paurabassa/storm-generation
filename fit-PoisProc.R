waves <- read.csv("./Clean.Data/waves-3hours.csv")
waves$Date <- as.POSIXct(strptime(waves$Date, "%Y-%m-%d %H:%M:%S"))
source("functions_pots.R")
pots.storm <- ExtractPOTPointProcess(waves$hs,2.5)

storm.start <- as.POSIXct(waves$Date[pots.storm$p.exc])
storm.end   <- as.POSIXct(waves$Date[pots.storm$p.exc + pots.storm$c.siz])

yearfraction <- function(date){
   print(date) 
   aux   <- as.POSIXlt(date)
   year  <- aux$year
   lyear <- difftime(as.POSIXct(strptime(paste(1901+year,"-01-01",sep=""),"%Y-%m-%d")),  
                     as.POSIXct(strptime(paste(1900+year,"-01-01",sep=""),"%Y-%m-%d")),
                     units="hours")
   fyear <- (24*aux$yday + aux$hour)/lyear[[1]]
}

ss.frac <- sapply(storm.start,yearfraction)
se.frac <- sapply(storm.end  ,yearfraction)
l <- length(ss.frac)
which(ss.frac[2:l]-se.frac[1:(l-1)]<0)

faux<- function(x){ if(x<0){ y <-x+1}else{ y<-x}; y}
hist(sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux),breaks=50)
nyears <- as.POSIXlt(tail(waves$Date,n=1))$year - as.POSIXlt(waves$Date[1])$year

lambda1 <- function(t, theta){
  theta[1] + 0*t
}

lambda2 <- function(t, theta) {
  theta[1] + theta[2] * sin(2*pi*t) + theta[3] * cos(2*pi*t)
} 

lambda3 <- function(t, theta) { 
  theta[1] + 
  theta[2] * dnorm(tan(pi*t+theta[3]), mean = theta[4], sd = theta[5])
}

x <- seq(-1,1,0.01)
y <- sapply(x,function(x) { lambda1(x,1)}) 
y <- sapply(x,function(x) { lambda2(x,c(1,0.5,0))})
y <- sapply(x,function(x) { lambda3(x,c(0,1,0,0,1))})

plot(x,y,type="l")

integrate(lambda1, 1, lower=0.1, upper=0.5)
integrate(lambda2, c(1,0.5,0), lower=0.1, upper=0.5)
integrate(lambda3, c(0,1,0,0,1), lower=0.1, upper=0.5)

my.integrate <- function(lambda, theta, lower, upper){
  if(upper-lower>0){
    I <- integrate(lambda, theta, lower=lower,upper=upper)[[1]]
  }else{
    I <- integrate(lambda, theta, lower = upper, upper = 1)[[1]] + 
         integrate(lambda, theta, lower = 0, upper = lower)[[1]]
  }
  I
}

ll.lam1 <- function(theta){ 
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


llike.lam2 <- function(theta1, theta2, theta3){ 
  theta <- c(theta1, theta2, theta3)
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

fit.lam3  <- mle( llike.lam3, method="BFGS", lower=c(0,0,0,-Inf,0),
                  start = list(theta1 = 1, theta2=1, theta3=0., 
                               theta4 = 0, theta5=1))
thet.lam3 <- coef(fit.lam3)

