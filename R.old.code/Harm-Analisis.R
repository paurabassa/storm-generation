set.seed(1234)

# read data 
p.clim <- read.csv("./Clean.Data/clim-hourly.csv")  # past climate
p.clim$Date <- as.POSIXct(strptime(p.clim$Date, "%Y-%m-%d %H:%M:%S","GMT"))
f.tide <- read.csv("./Clean.Data/future-tides.csv") # future tide
f.tide$Date <- as.POSIXct(strptime(f.tide$Date, "%Y-%m-%d %H:%M:%S", "GMT"))


# set variable of study in a logaritmic scale 
log.hs <- log(p.clim$hs) 
# load code for Harmonic Analysis
source("R.functions/Harmonic.R")

# split signal in log.hs into periodic and residual components 
split.lg.hs <- SubstractPeriodSignal(log.hs, c(1/(365.25*24),1/(8)), c(15,5))
ts.orig <- log.hs  #original time series
ts.peri <- split.lg.hs$periodic
ts.stat <- split.lg.hs$residual
p.clim$lg.hs.p <- split.lg.hs$periodic
p.clim$lg.hs.r  <- split.lg.hs$residual

##############################################################
# ------ study of the point process of excedances -----------#
##############################################################


# ---- compute the threshold of the stationary -----------------
# ------ part to compute all storm candidates ------------------
thres.orig <- log(2.5)                          #original threshold
thres.stat <- thres.orig - max(ts.peri)         #candidates threshold

source("R.functions/ExtractPOTpp.R")
pots <- ExtractPOTPointProcess(ts.stat,thres.stat)

storms <- p.clim[ p.clim$lg.hs.r > thres.stat,]
str(pots)
true.exc <-storms$lg.hs.p+storms$lg.hs.r > thres.orig
plot(storms$lg.hs.p, storms$lg.hs.r, col=ifelse(true.exc,"red","black"))
lines(seq(-1.3,0.1,0.01),thres.orig -seq(-1.3,0.1,0.01))

#############################################################################
# -- computation of the empirical distribution function with smoothing ---- #
# ------------------ for the storms inter arrival times ------------------- #
#  #
############################################################################

sample <- pots$nex.exc
dens.lg   <- density(log(sample))
plot(dens.lg)

xfit <- exp(dens.lg$x)
yfit <- dens.lg$y*exp(-dens.lg$x)

#load auxirial function CDF (cumulative distribution function)
source("R.functions/funcs_HA.R")
#function to sample the new distrib
r.Approx.Empirical <- function(nsample=1, xfit,yfit){
  approx(CDF(xfit,yfit),xfit,runif(nsample))[[2]]
}

resample <- r.Approx.Empirical(length(sample), xfit,yfit)
new.nex.exc <- as.integer(resample-min(resample)+1)


op <- par(mfrow=c(2,2))
h    <- hist(sample,probability=T,breaks=100)
#lines(xfit,yfit,xlim=c(0,max(sample)),col="red")
plot(h$mids,h$density,log='xy')
lines(xfit,yfit, col="red",type='l')

h <- hist(new.nex.exc,breaks=100,prob=T)
plot(h$mids, h$density, log='xy')
lines(xfit,yfit, col="red",type='l')
par(op)



####################################################
# Part 3: Generate a sample                        #
####################################################
FindStorm <- function(nex.exc, pots){
  aux <- which(pots$nex.exc == nex.exc)
  if (length(aux)==0){
    i.aux <-which.min(abs(pots$nex.exc - nex.exc))
    aux <- which(pots$nex.exc== pots$nex.exc[i.aux])
    index <- aux[sample(1:length(aux),1)]
  } else {
    index <- aux[sample(1:length(aux),1)]
  }
  index
}

# storms generation
n.str  <- as.POSIXct("20150101 00:00",format ="%Y%m%d %H:%M") #next strom time
i      <- 1
future <- as.POSIXct("20650101 00:00",format ="%Y%m%d %H:%M")
s.t   <- vector()  # storm start time
s.ind <- vector()  # storm index to be cloned
gap   <- vector()  # time to the next storm

while(n.str < future){
  s.t[i]   <- as.POSIXct(n.str, origin="1970-01-01 00:00.00 UTC")
#  yf       <- yearfraction(date=n.str)         #year fraction
#  n.gap.aux<- r.pp.cond(lambda3,thet.lam3,yf) #gap in year units
  gap[i]   <- max(as.integer(r.Approx.Empirical(1, xfit,yfit)),1)
  s.ind[i] <- FindStorm(gap[i],pots)               #storm index
  n.dur    <- pots$c.siz[s.ind[i]]  #storm duration hours
  n.str    <- as.POSIXct(s.t[i] + 3600*(gap[i] + n.dur), 
                         origin="1970-01-01 00:00.00 UTC") #time to next storm
  i        <- i+1
}

new.storms <- data.frame(s.times= as.POSIXct(s.t,
                                          origin="1970-01-01 00:00.00 UTC"),
                         s.index= s.ind, s.dur= pots$c.siz[s.ind],
                         s.gap=gap)
length(pots[[1]])/11
length(new.storms[[1]])/50

storm.climate <- data.frame(Date = as.POSIXct(vector(mode="numeric"),
                                           origin= '1970-01-01 00:00.00 UTC'),
                        hs = vector(mode="numeric"), fp  = vector(mode="numeric"),
                        tm = vector(mode="numeric"), dir = vector(mode="numeric"),
                        U10 = vector(mode="numeric"), V10 = vector(mode="numeric"),
                        a.tide=vector(mode="numeric"), res= vector(mode="numeric"))
i.aux<-1
for(i in 1:length(new.storms$s.times)){
  print(paste(i," out of ",length(new.storms$s.times),sep=""))
  for(j in 1:pots$c.siz[new.storms$s.ind[i]]){
    i.past   <- pots$p.exc[new.storms$s.ind[i]] + j-1
    date.fut <- as.POSIXct(new.storms$s.times[i] + 3600*(j-1),
                                          origin="1970-01-01 00:00.00 UTC")
#    index.fut.tide <- which(f.tide$Date == date.fut

    tide.fut <- f.tide[which(f.tide$Date == date.fut),]$Level
    aux.date <- paste("2004-",format(date.fut, "%m-%d %H:%M"), sep="")
    date.hs.p <- as.POSIXct(strptime(aux.date, "%Y-%m-%d %H:%M", "GMT"))
    lg.hs.p <- p.clim[which(p.clim$Date == date.hs.p),"lg.hs.p"]
    lg.hs.r <- p.clim[i.past,"lg.hs.r"]
    storm.climate[i.aux,"hs"] = exp(lg.hs.p + lg.hs.r)
    storm.climate[i.aux,"Date"]    = date.fut
    storm.climate[i.aux, c("fp", "tm", "dir", "U10", "V10", "res")] =
                          p.clim[i.past, c("fp", "tm", "dir", "U10", "V10", "res")]
    storm.climate[i.aux, "a.tide"] = tide.fut
    i.aux <- i.aux + 1
  }
}

true.storms <- storm.climate[storm.climate$hs>2.5,]
write.csv(true.storms, "New.Data/storms-Harm.csv",row.names=F)

orig.storms <- p.clim[p.clim$hs>2.5,]
write.csv(orig.storms, "New.Data/storms-Orig.csv",row.names=F)






################################################
#
# IGNORE FROM HERE ONWARDS
#
##########################















#####################################################################
# ---- select stroms with same inter.arr.times as the ones 
# ---- saved in resample
###################################################


FindStorm(  1, pots)
FindStorm(320, pots)
pots$nex.exc[FindStorm(1,pots)]

# generate a new list of Peaks over threshold, in which the inter arrivals 
# times are the ones of resample, p.exc contains the start of the mirror 
# storm in data.in, and c.siz contains the duration of the mirror storm 

nresample <- length(resample)

new.pots <- list( p.exc = vector(mode="numeric", length = nresample),
                  c.siz = vector(mode="numeric", length = nresample),
                  nex.exc = new.nex.exc)
#set.seed(1234)
for(i in 1:nresample){
  i.pots.storm <-  FindStorm(new.pots$nex.exc[i],pots)
  new.pots$p.exc[i] = pots$p.exc[i.pots.storm]
  new.pots$c.siz[i] = pots$c.siz[i.pots.storm]
}

plot(new.pots$c.siz,new.pots$p.exc)

############################################################
#  auxiliary function to print the storms out of a climate #
############################################################

ExtractClimateStorms <- function(climate){
  storm.climate <- list(Date = as.POSIXct(vector(mode="numeric"),
                                           origin= '1970-01-01 00:00.00 UTC'),
                        hs = vector(mode="numeric"), fp  = vector(mode="numeric"),
                        tm = vector(mode="numeric"), dir = vector(mode="numeric"),
                        U10 = vector(mode="numeric"), V10 = vector(mode="numeric"))
  i.aux<-1
  for(i in 1:length(climate$hs)){
#  for(i in 1:3000){ 
    if(climate$hs[i]>2.5){
      storm.climate$Date[i.aux] = climate$Date[i]
      storm.climate$hs[i.aux]   = climate$hs[i]
      storm.climate$fp[i.aux]   = climate$fp[i]
      storm.climate$tm[i.aux]   = climate$tm[i]
      storm.climate$dir[i.aux]  = climate$dir[i]
      storm.climate$U10[i.aux]  = climate$U10[i]
      storm.climate$V10[i.aux]  = climate$V10[i]
      i.aux <- i.aux + 1
    }
  }
  as.data.frame(storm.climate)
}

storm.orig <- ExtractClimateStorms(p.clim)
str(storm.orig)
plot(storm.orig$Date, storm.orig$hs)
write.csv(storm.orig, "orig-storms.csv",row.names=F)


########################################################
#  generate data like clim.in for the resampled storms #
########################################################

nstorms <- nresample        # total numbe of resampled storms
#nstorms <- 1000
max.time <- length(ts.peri)
max.time <- 8*365.25*20          # 4 years
set.seed(1234)

new.clim <- data.frame(Date = as.POSIXct(vector(mode="numeric",length=max.time),
                                           origin= '1970-01-01 00:00.00 UTC'),
                 hs   = vector(mode="numeric", length=max.time),
                 fp   = vector(mode="numeric", length=max.time),
                 tm   = vector(mode="numeric", length=max.time),
                 dir  = vector(mode="numeric", length=max.time),
                 U10  = vector(mode="numeric", length=max.time),
                 V10  = vector(mode="numeric", length=max.time))
new.clim$Date[1] <- as.POSIXct("2013-01-01 GMT")
for(i in 2:max.time) new.clim$Date[i]  <-  new.clim$Date[i-1] + 1*3600

i.time <- 1
i.storm <- 1
while(i.time +  new.pots$c.siz[i.storm] < max.time & i.storm < nstorms){
#  start.storm <- new.pots$p.exc[i.storm]
#  duration    <- new.pots$c.siz[i.storm]
  storm.dur   <- new.pots$c.siz[i.storm]
  storm.times <- new.pots$p.exc[i.storm]:(new.pots$p.exc[i.storm] +
                                          new.pots$c.siz[i.storm])
  new.clim$hs[i.time:(i.time+storm.dur)] <-
                exp(ts.peri[i.time:(i.time+storm.dur)] + ts.stat[storm.times])
  new.clim$fp[i.time:(i.time+storm.dur)] <- p.clim$fp[storm.times]
  new.clim$tm[i.time:(i.time+storm.dur)] <- p.clim$tm[storm.times]
  new.clim$dir[i.time:(i.time+storm.dur)] <- p.clim$dir[storm.times]
  new.clim$U10[i.time:(i.time+storm.dur)] <- p.clim$U10[storm.times]
  new.clim$V10[i.time:(i.time+storm.dur)] <- p.clim$V10[storm.times]

  i.time  <- i.time + storm.dur + new.pots$nex.exc[i.storm]
  i.storm <- i.storm + 1
}
                                                                         
i.time
i.storm

new.storms <- ExtractClimateStorms(new.clim)
str(new.storms)
plot(new.storms$Date, new.storms$hs)
write.csv(new.storms, "storms-new.csv",row.names=F)

