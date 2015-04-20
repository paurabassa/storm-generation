#R file to analyise the clustering of storm events. 

##################################################
#  Part 1: Preliminary analysis of the Hs signal #
##################################################
#read ocean climate data (waves + tides)
clim <- read.csv("./Data.Processed/clim-hourly.csv")
#fix the format of the Date variable
clim$Date <- as.POSIXct(strptime(clim$Date, "%Y-%m-%d %H:%M:%S"))
str(clim)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0),bg="white")
plot(clim$Date, clim$hs, type='l', main="Original Time Series",
      xlab="Year", ylab="Significant Wave Height")
abline(a=2.5,b=0, col="red")
#pl1 <-recordPlot()

fil.1m <- rep(1,24*30)/(24*30)
fil.3m <- rep(1,24*30*3)/(24*30*3)
plot(clim$Date, log(clim$hs), type='l', main="Transformed Time Series", 
     xlab="Year", ylab="log(Height)", ylim=c(-4,3.5))
lines(clim$Date, filter(log(clim$hs), fil.1m), col="red")
lines(clim$Date, filter(log(clim$hs), fil.3m), col="blue")
legend("topleft",
       c("observed", "one month average", "three month average"), 
       lty=1, col=c("black","red","blue"))
#pl2 <-recordPlot()

lg.hs.acf <- acf(log(clim$hs),lag.max=365.25*4*24, xlab="Lag (hours)", 
                 main="Auto Correlation Function")
#pl3 <- recordPlot()
 
spectrum(log(clim$hs),main="Raw Periodogram")

dev.print(device=png, file="Figures/1-hs.png",width=800,height=800)

#pl4 <- recordPlot()


###############################################
#    Part 2: Clusering of storms              #
###############################################


#select the climate of stroms 
clim.storms <- clim[clim$hs>2.5,]
plot(clim.storms[,c("hs", "fp", "tm", "dir", "res")])

#preliminary plot of the behaviour of clusters
storms.times <- as.POSIXlt(clim.storms$Date)
plot(storms.times$yday,storms.times$year+1900,ylim=c(2013,1990))

#
rng <- range(storms.times$yday)
inc <- (rng[2]-rng[1])/12
seq(rng[1],rng[2],inc)
hist(storms.times$yday,breaks= seq(rng[1],rng[2],inc))

# to use extended time lapse (where tides are not represented) we can do:

waves <- read.csv("./Data.Processed/waves-3hours.csv")
waves$Date <- as.POSIXct(strptime(waves$Date, "%Y-%m-%d %H:%M:%S"))
clim.storms <- waves[waves$hs>2.5,]
storms.times <- as.POSIXlt(clim.storms$Date)
plot(storms.times$yday,storms.times$year+1900,ylim=c(2013,1980))
rng <- range(storms.times$yday)
inc <- (rng[2]-rng[1])/12
seq(rng[1],rng[2],inc)
hist(storms.times$yday,breaks= seq(rng[1],rng[2],inc))

faux <- function(x){
  if(x<183){ y<- x+183
  }else{     y<- x-183}
  y
}
aux <- sapply(storms.times$yday, faux)
hist(aux,breaks= seq(rng[1],rng[2],inc))


# Storm times using Peaks over threshold and clusteriing 

source("R.functions/ExtractPOTpp.R")
pots.storm <- ExtractPOTPointProcess(waves$hs,2.5)
storms.times <- as.POSIXlt(waves$Date[pots.storm$p.exc])
storms.durat <- pots.storm$c.siz

sct <- data.frame(Day=storms.times$yday, 
                           Year=storms.times$year+1900, 
                           Duration=3*storms.durat)
symbols(x=sct$Day, y=sct$Year, circles=sct$Duration, inches=1/9, ann=F, bg="steelblue2", fg=NULL)

library("ggplot2")
qplot(x=Day, y=Year, data=sct, size=Duration)
