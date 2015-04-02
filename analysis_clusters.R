#R file to analyise the clustering of storm events. 

#read ocean climate data (waves + tides)
#clim <- read.csv("./Clean.Data/clim-3hours.csv")
clim <- read.csv("./Clean.Data/clim-hourly.csv")

#fix the format of the Date variable
clim$Date <- as.POSIXct(strptime(clim$Date, "%Y-%m-%d %H:%M:%S"))
str(clim)

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

waves <- read.csv("./Clean.Data/waves-3hours.csv")
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
symbols(x=sct$ev1, y=sct$ev2, circles=sct$ev3, inches=1/9, ann=F, bg="steelblue2", fg=NULL)

library("ggplot2")
qplot(x=Day, y=Year, data=sct, size=Duration)
