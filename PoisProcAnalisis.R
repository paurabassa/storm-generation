set.seed(1234)
##
##
##
##
##
##
##

####################
# Part 1: Data fit #
####################

# Read data
waves <- read.csv("./Clean.Data/waves-hourly.csv")
waves$Date <- as.POSIXct(strptime(waves$Date, "%Y-%m-%d %H:%M:%S","GMT"))

# Extract the point process of excedances (storms)
source("R.functions/ExtractPOTpp.R")
pots.storm <- ExtractPOTPointProcess(waves$hs,2.5)
storm.start <- as.POSIXct(waves$Date[pots.storm$p.exc])
storm.end   <- as.POSIXct(waves$Date[pots.storm$p.exc + pots.storm$c.siz])

# Computation of the year fraction times and sanity test
source("R.functions/YearFraction.R")
ss.frac <- sapply(storm.start, yearfraction)
se.frac <- sapply(storm.end  , yearfraction)
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
plot(sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux),pots.storm$c.siz[1:(l-1)])



# maximum likelihood fits of the gaps to a pp.distribution with intensity 
# function lambda_i, i=1,2,3. For lambda_3 it is necessary to estimate 3
# parameters first, and then add the other 2. 
source("R.functions/funcs_PPA.R")

fit.lam1  <- mle(ll.lam1,list(theta=1))
thet.lam1 <- coef(fit.lam1)[[1]]

fit.lam2  <- mle( llike.lam2, start = list(theta1 = 1, theta2=0.4, theta3=0.4))
thet.lam2 <- coef(fit.lam2)

fit.lam3.aux <- mle( llike.lam3.aux, method="L-BFGS-B", lower=c(0.01,0.01,0),
                  start = list(theta1 = 0, theta2=1, theta3=0.))
thet.lam3.aux <- coef(fit.lam3.aux)
fit.lam3  <- mle( llike.lam3, method="BFGS", lower=c(0,0.01,0,-Inf,0),
                  start = list(theta1 = thet.lam3.aux[[1]], theta2=thet.lam3.aux[[2]], 
                               theta3 = thet.lam3.aux[[3]], theta4 = 0, theta5=1))
thet.lam3 <- coef(fit.lam3)

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

hist(ss.frac, prob = T)
lines(x, y1, lty=1)
lines(x, y2, lty=2)
lines(x, y3, lty=3)

# quantile-quantile plots of the empirical inter arrival times vs 
# fitted Poission Process

par(mfrow=c(1,3))

x.lam1 <- vector("numeric")
for (i in 1:(l-1)) x.lam1[i]<- p.pp.cond(gaps[i],lambda1, thet.lam1, se.frac[i])
qqplot(qunif(ppoints(500)), x.lam1)
lines(seq(-0.1,1.1,0.0001),seq(-.1,1.1,0.0001))

# lambda2
x.lam2 <- vector("numeric")
for (i in 1:(l-1)) x.lam2[i]<- p.pp.cond(gaps[i],lambda2, thet.lam2, se.frac[i])
qqplot(qunif(ppoints(500)), x.lam2)
lines(seq(-0.1,1.1,0.0001),seq(-.1,1.1,0.0001))

# lambda3 
x.lam3 <- vector("numeric")
for (i in 1:(l-1)) x.lam3[i]<- p.pp.cond(gaps[i],lambda3, thet.lam3, se.frac[i])
qqplot(qunif(ppoints(500)), x.lam3)
lines(seq(-0.1,1.1,0.0001),seq(-.1,1.1,0.0001))


#############################
# Part 3: Generate a sample #
#############################

# Given a set of indexed storms from 1:nst we generate a new resample of 
# storms, which is represented by
# storm start
# storm index 
# storm duration

# storms generation
n.str  <- as.POSIXct("20150101 00:00",format ="%Y%m%d %H:%M") #next strom time
i      <- 1   
future <- as.POSIXct("20650101 00:00",format ="%Y%m%d %H:%M")
n.st   <- vector()
n.sind <- vector()
n.gap  <- vector()

while(n.str < future){ 
  n.st[i]  <- as.POSIXct(n.str, origin="1970-01-01 00:00.00 UTC") 
  yf       <- yearfraction(date=n.str)         #year fraction
  n.gap.aux<- r.pp.cond(lambda3,thet.lam3,yf) #gap in year units
  n.gap[i] <- as.integer(n.gap.aux*365*24) #gap in seconds truncated at hours
  n.sind[i]<- sample(1:l,size=1)               #storm index
  n.dur    <- pots.storm$c.siz[n.sind[i]]  #storm duration in seconds 
  n.str    <- as.POSIXct(n.st[i] + 3600*(n.gap[i] + n.dur), origin="1970-01-01 00:00.00 UTC")
  i        <- i+1
}

# save new set of stroms into a data set
new.storms <- data.frame(s.times= as.POSIXct(n.st, 
                                          origin="1970-01-01 00:00.00 UTC"), 
                         s.index= n.sind, s.dur= pots.storm$c.siz[n.sind], 
                         s.gap=n.gap)

# qqplot of generated data vs original set of storms
qqplot(new.storms$s.gap, as.integer(365*24*gaps))
lines(seq(-10,10000,0.1), seq(-10,10000,0.1))

#careful!: 
#storms x year original data:
length(pots.storm[[1]])/nyears
 
#storms x year simulated: 
length(new.storms[[1]])/50


## 
##  Write the resample of future storms with its corresponding predicted tide.
## 

p.clim <- read.csv("./Clean.Data/clim-hourly.csv")  # past climate
p.clim$Date <- as.POSIXct(strptime(p.clim$Date, "%Y-%m-%d %H:%M:%S","GMT"))
f.tide <- read.csv("./Clean.Data/future-tides.csv") # future tide
f.tide$Date <- as.POSIXct(strptime(f.tide$Date, "%Y-%m-%d %H:%M:%S", "GMT"))


storm.climate <- data.frame(Date = as.POSIXct(vector(mode="numeric"),
                                           origin= '1970-01-01 00:00.00 UTC'),
                        hs = vector(mode="numeric"), fp  = vector(mode="numeric"),
                        tm = vector(mode="numeric"), dir = vector(mode="numeric"),
                        U10 = vector(mode="numeric"), V10 = vector(mode="numeric"),
                        a.tide=vector(mode="numeric"), res= vector(mode="numeric"))
i.aux<-1
for(i in 1:length(new.storms$s.times)){
  print(paste(i," out of ",length(new.storms$s.times),sep=""))
  for(j in 1:pots.storm$c.siz[new.storms$s.ind[i]]){
    i.past   <- pots.storm$p.exc[new.storms$s.ind[i]] + j-1
    date.fut <- as.POSIXct(new.storms$s.times[i] + 3600*(j-1),
                                          origin="1970-01-01 00:00.00 UTC")
#    index.fut.tide <- which(f.tide$Date == date.fut
    
    tide.fut <- f.tide[which(f.tide$Date == date.fut),]$Level
      
    storm.climate[i.aux,"Date"]    = date.fut
    storm.climate[i.aux, c("hs", "fp", "tm", "dir", "U10", "V10", "res")] =
                          p.clim[i.past, c("hs", "fp", "tm", "dir", "U10", "V10", "res")]
    storm.climate[i.aux, "a.tide"] = tide.fut
    i.aux <- i.aux + 1
  }
}


write.csv(storm.climate, "New.Data/storms-Pois.csv",row.names=F)

