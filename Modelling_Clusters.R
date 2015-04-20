##
## Modelling of clustering of storms through different approxes. 
## The main random variable to analyse is the interarrival time 
## of storms. First we address it via homogenous Poisson Process, a 
## non-homogenuous one, and via a GAM model, using the year fraction 
## as explanatory variable. Also includes an approximation of the 
## density function via gausian kernels after the periodic variable is 
## removed. 
##

#######################
# Part 0 : preliminary checks and data preparation
#######################

#set random seed
set.seed(1234)

# Read data
clim       <- read.csv("./Data.Processed/clim-hourly-ext.csv")
clim$Date  <- as.POSIXct(strptime(waves$Date, "%Y-%m-%d %H:%M:%S","GMT"))
waves      <- read.csv("./Data.Processed/waves-hourly.csv")
waves$Date <- as.POSIXct(strptime(waves$Date, "%Y-%m-%d %H:%M:%S","GMT"))

# Extract the point process of excedances (storms)
source("R.functions/ExtractPOTpp.R")
pots.storm <- ExtractPOTPointProcess(waves$hs,2.5)
storm.start <- as.POSIXct(waves$Date[pots.storm$p.exc])
storm.end   <- as.POSIXct(waves$Date[pots.storm$p.exc + pots.storm$c.siz])
summary(difftime(storm.end,storm.start))

# Computation of the year fraction times and sanity test
source("R.functions/YearFraction.R")
ss.frac <- sapply(storm.start, yearfraction)
se.frac <- sapply(storm.end  , yearfraction)
l <- length(ss.frac)
which(ss.frac[2:l]-se.frac[1:(l-1)]<0)

# histogram of diferences taking into acount the jump from one year to the next
faux<- function(x){ if(x<0){ y <-x+1}else{ y<-x}; y}
gaps <- sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux)
hist(gaps,breaks=25)
nyears <- as.POSIXlt(tail(waves$Date,n=1))$year - as.POSIXlt(waves$Date[1])$year +1

# check if there is a correlation between storm duration (c.siz) and the
# inter arrival times
cor(sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux),pots.storm$c.siz[1:(l-1)],method="k")
plot(rank(sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux),ties.method="r")/l,rank(pots.storm$c.siz[1:(l-1)],ties.method="r")/l)

# check if there is a correlation between the year fraction and the 
# inter arrival times
cor(sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux),ss.frac[2:l],method="k")
plot(rank(sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux),ties.method="r")/l,rank(ss.frac[2:l],ties.method="r")/l)

# manual fit to weibull and exponential distributions, just for check
plot(sort(gaps/mean(gaps)),1-ppoints(l-1,0),log="y")
lines(sort(gaps/mean(gaps)), exp(-sort(gaps/mean(gaps))))
lines(sort(gaps/mean(gaps)), 1-pweibull(sort(gaps/mean(gaps)),.55))

####################
# Part 1: Data fit #
####################

# maximum likelihood fits of the gaps to a pp.distribution with intensity
# function lambda_i, i=1,2,3. For lambda_3 it is necessary to estimate 3
# parameters first, and then add the other 2.
source("R.functions/funcs_PPA.R")

# homogenous Poisson
fit.lam1  <- mle(ll.lam1,list(theta=1))
thet.lam1 <- coef(fit.lam1)[[1]]

# non homogenous Poisson 1
fit.lam2  <- mle( llike.lam2, start = list(theta1 = 1, theta2=0.4, theta3=0.4))
thet.lam2 <- coef(fit.lam2)

# non homogenous Poisson 2
fit.lam3.aux <- mle( llike.lam3.aux, method="L-BFGS-B", lower=c(0.01,0.01,0),
                  start = list(theta1 = 0, theta2=1, theta3=0.))
thet.lam3.aux <- coef(fit.lam3.aux)
fit.lam3  <- mle( llike.lam3, method="L-BFGS-B", lower=c(0,0.01,0,-Inf,0),
                  start = list(theta1 = thet.lam3.aux[[1]], theta2=thet.lam3.aux[[2]],
                               theta3 = thet.lam3.aux[[3]], theta4 = 0, theta5=1))
thet.lam3 <- coef(fit.lam3)

# GAM modelling 
library(gamlss)
tt <- ss.frac[1:(l-1)]

fit.1 <- gamlss(gaps   ~  sin(2*pi*tt) + cos(2*pi*tt) + sin(4*pi*tt) + cos(4*pi*tt) 
#                          + sin(8*pi*tt) + cos(8*pi*tt)
     , sigma.formula = ~  sin(2*pi*tt) + cos(2*pi*tt) + sin(4*pi*tt) + cos(4*pi*tt) 
#                          + sin(8*pi*tt) + cos(8*pi*tt)
                , family="BE")
centiles(fit.1,tt ,ylim=c(0,1))
wp(fit.1)
plot(fit.1)

# ??? kernel approximation to the  empirical distribution ??? 


############################
# Part 2: Diagnostic Plots #
############################

# Densities of the 3 fitted distributions vs empirical histogram
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

hist(se.frac, prob = T)
lines(x, y1, lty=1)
lines(x, y2, lty=2)
lines(x, y3, lty=3)

# probability-probability plots of the empirical inter arrival times vs
# fitted models

op <- par(mfrow=c(2,2),mar=c(4,4,1,1))

# lambda2
x.lam1 <- vector("numeric")
for (i in 1:(l-1)) x.lam1[i]<- p.pp.cond(gaps[i],lambda1, thet.lam1, se.frac[i])
qqplot(ppoints(500), x.lam1)
abline(0,1)

# lambda2
x.lam2 <- vector("numeric")
for (i in 1:(l-1)) x.lam2[i]<- p.pp.cond(gaps[i],lambda2, thet.lam2, se.frac[i])
qqplot(ppoints(500), x.lam2)
abline(0,1)

# lambda3
x.lam3 <- vector("numeric")
for (i in 1:(l-1)) x.lam3[i]<- p.pp.cond(gaps[i],lambda3, thet.lam3, se.frac[i])
qqplot(ppoints(500), x.lam3)
abline(0,1)

# GAM
para <- predictAll(fit.1)
p.BE <- pBE(gaps,para$mu, para$sigma)
qqplot(ppoints(500), p.BE)
abline(0,1)

par(op)


