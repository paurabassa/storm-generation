#R file to analyise the clustering of storm events. 

##################################################
#  Part 1: Preliminary analysis of the Hs signal #
##################################################
#read ocean climate data (waves + tides)
clim <- read.csv("./Data.Processed/clim-hourly.csv")
#fix the format of the Date variable
clim$Date <- as.POSIXct(strptime(clim$Date, "%Y-%m-%d %H:%M:%S"))
str(clim)

#plot original time series
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0),bg="white")
plot(clim$Date, clim$hs, type='l', main="Original Time Series",
      xlab="Year", ylab="Significant Wave Height")
abline(a=2.5,b=0, col="red")

#plot log of TS with rolling averages
fil.1m <- rep(1,24*30)/(24*30)
fil.3m <- rep(1,24*30*3)/(24*30*3)
plot(clim$Date, log(clim$hs), type='l', main="Transformed Time Series", 
     xlab="Year", ylab="log(Height)", ylim=c(-4,3.5))
lines(clim$Date, filter(log(clim$hs), fil.1m), col="red")
lines(clim$Date, filter(log(clim$hs), fil.3m), col="blue")
legend("topleft",
       c("observed", "one month average", "three month average"), 
       lty=1, col=c("black","red","blue"))

#plot acf
lg.hs.acf <- acf(log(clim$hs),lag.max=365.25*4*24, xlab="Lag (hours)", 
                 main="Auto Correlation Function")
#plot spectrum
spectrum(log(clim$hs),main="Raw Periodogram")

#save file
dev.print(device=png, file="Figures/1-hs.png",width=800,height=800)



###############################################
#    Part 2: Clusering of storms              #
###############################################
library("ggplot2")
source("R.functions/ExtractPOTpp.R")

pots.storm   <- ExtractPOTPointProcess(clim$hs,2.5)
storms.times <- as.POSIXlt(clim$Date[pots.storm$p.exc])
storms.dur   <- pots.storm$c.siz

sct <- data.frame(Day=storms.times$yday,
                  Year=as.factor(storms.times$year+1900),
                  Duration=storms.dur)
require(gridExtra)

Month <- as.factor(storms.times$mon)
levels(Month) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                "Sep", "Oct", "Nov", "Dec")

p1 <- qplot(Month)
p2 <- qplot(x=Day, y=Year, data=sct, size=Duration)
grid.arrange(p1,p2,ncol=2)

dev.print(device=png, file="Figures/2-clusters.png",width=800,height=400)

###############################################################
# Part 3: Same analysis after removing periodic variable      #
###############################################################
source("R.functions/Harmonic.R")
str(clim)
log.hs       <- log(clim$hs)
split.lg.hs  <- SubstractPeriodSignal(log.hs, c(1/(365.25*24),1/(8)), c(15,5))
clim$lg.hs.p <- split.lg.hs$periodic
clim$lg.hs.r <- split.lg.hs$residual
write.csv(clim,"Data.Processed/clim-hourly-ext.csv")

# compute the threshold of the stationary part
# to compute all storm candidates
thres.orig <- log(2.5)                          #original threshold
thres.stat <- thres.orig - max(clim$lg.hs.p)  #candidates threshold

pots.storm <- ExtractPOTPointProcess(clim$lg.hs.r, thres.stat)

storms.times <- as.POSIXlt(clim$Date[pots.storm$p.exc])
storms.dur   <- pots.storm$c.siz

sct <- data.frame(Day=storms.times$yday,
                  Year=as.factor(storms.times$year+1900),
                  Duration=storms.dur)

Month <- as.factor(storms.times$mon)
levels(Month) <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
                "Sep", "Oct", "Nov", "Dec")

p1 <- qplot(Month)
p2 <- qplot(x=Day, y=Year, data=sct, size=Duration)
grid.arrange(p1,p2,ncol=2)

dev.print(device=png, file="Figures/3-clusters-B.png",width=800,height=400)

#below here plot dependency
storms <- clim[ clim$lg.hs.r > thres.stat,]
str(pots)
true.exc <-storms$lg.hs.p+storms$lg.hs.r > thres.orig
plot(storms$lg.hs.p, storms$lg.hs.r, col=ifelse(true.exc,"red","black"))
lines(seq(-1.3,0.1,0.01), thres.orig -seq(-1.3,0.1,0.01))
dev.print(device=png, file="Figures/4-exedance-harm.png",width=600,height=600)

