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

interv <- range(p.clim$lg.hs.p)
wid <- abs(interv[1]-interv[2])/100

# storms generation
n.str  <- as.POSIXct("20150101 00:00",format ="%Y%m%d %H:%M") #next strom time
i      <- 1
future <- as.POSIXct("20650101 00:00",format ="%Y%m%d %H:%M")
s.t   <- vector()  # storm start time
s.ind <- vector()  # storm index to be cloned
gap   <- vector()  # time to the next storm

while(n.str < future){
  s.t[i]    <- as.POSIXct(n.str, origin="1970-01-01 00:00.00 UTC") #storm date
  gap[i]   <- max(as.integer(r.Approx.Empirical(1, xfit,yfit)),1) 
  #time to next storm
  
# select a storm with similar value of  lg.hs.p for the corresponding
  aux.date <- paste("2004-",format(n.str, "%m-%d %H:%M"), sep="")
  date.hs.p <- as.POSIXct(strptime(aux.date, "%Y-%m-%d %H:%M", "GMT"))
  lg.hs.p <- p.clim[which(p.clim$Date == date.hs.p),"lg.hs.p"] #value of lg.hs.p at the start of the storm 
  s.ind[i] <- sample(which(abs(p.clim$lg.hs.p[pots$p.exc] - lg.hs.p)<wid),1)  
#  yf       <- yearfraction(date=n.str)         #year fraction
#  n.gap.aux<- r.pp.cond(lambda3,thet.lam3,yf) #gap in year units
#  s.ind[i] <- FindStorm(gap[i],pots)               #storm index
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
write.csv(true.storms, "New.Data/storms-Harm-2.csv",row.names=F)


