##
##  Script to generate n.s samples of storms under different 
##  climate change scenarios rcp26 and rpc85. Also generates 
##  new samples with different cluster scenarios. 
##
##
##

set.seed(123456)
n.s = 3 # this is the number of samples for each scenario

#############################################################
#     Main auxiliary function to generate storm samples     #
#############################################################

GenerateStormSample <- function(rGap, rSind, pots.storm, p.clim, f.tide){ 
  #function that generates a resample of storms
  # rGap       -> input function that generates the gap to the next storm
  # rSind      -> input function that generates the storm index
  # pots.storm -> input dataframe with a catalog of storms
  # p.clim     -> input dataframe with the past climate
  # f.tide     -> input dataframe with the future tides
  n.str  <- as.POSIXct("20150101 00:00",format ="%Y%m%d %H:%M") #next strom time
  i      <- 1
  future <- as.POSIXct("20650101 00:00",format ="%Y%m%d %H:%M")
  n.st   <- vector()
  n.sind <- vector()
  n.gap  <- vector()

  while(n.str < future){
    n.st[i]  <- as.POSIXct(n.str, origin="1970-01-01 00:00.00 UTC")
    yf      <<- yearfraction(date=n.str)         #year fraction
    n.gap[i] <- rGap() 
    n.sind[i]<- rSind()           #storm index
    n.dur    <- pots.storm$c.siz[n.sind[i]]  #storm duration in seconds
    n.str    <- as.POSIXct(n.st[i] + 3600*(n.gap[i] + n.dur), 
                           origin="1970-01-01 00:00.00 UTC")
    i        <- i+1
  }

  # save new set of stroms into a data set
  # this data frame is unnecessary but makes the code more readable
  new.storms <- data.frame(s.times= as.POSIXct(n.st, 
                                        origin="1970-01-01 00:00.00 UTC"),
                         s.index= n.sind, s.dur= pots.storm$c.siz[n.sind],
                         s.gap=n.gap) 
  
  # this is the data frame that we will return 
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
      if(date.fut >=future){ break} 
  # this line below slows down the code a lot
      tide.fut <- f.tide[which(f.tide$Date == date.fut),]$Level

      storm.climate[i.aux,"Date"]    = date.fut
      storm.climate[i.aux, c("hs", "fp", "tm", "dir", "U10", "V10", "res")] =
           p.clim[i.past, c("hs", "fp", "tm", "dir", "U10", "V10", "res")]
      storm.climate[i.aux, "a.tide"] = tide.fut
      i.aux <- i.aux + 1
    }
  }
  nstorms <<- length(new.storms$s.times)
  storm.climate
} 


#####################################
#  Load data and 1st two scenarios  #
#  as a test to the function above  #
#####################################

# Read data
clim       <- read.csv("./Data.Processed/clim-hourly-ext.csv")
clim$Date  <- as.POSIXct(strptime(clim$Date, "%Y-%m-%d %H:%M:%S","GMT"))
f.tide <- read.csv("./Data.Processed/future-tides.csv")
f.tide$Date <-  as.POSIXct(strptime(f.tide$Date,"%Y-%m-%d %H:%M:%S","GMT"))


# Extract Storm
source("R.functions/ExtractPOTpp.R")
source("R.functions/YearFraction.R")
pots.storm <- ExtractPOTPointProcess(clim$hs,2.5)
storm.start <- as.POSIXct(clim$Date[pots.storm$p.exc])
storm.end   <- as.POSIXct(clim$Date[pots.storm$p.exc + pots.storm$c.siz])
ss.frac <- sapply(storm.start, yearfraction)
se.frac <- sapply(storm.end  , yearfraction)

yf      <-0
nstorms <- 0

l <- length(ss.frac)
nyears <- as.POSIXlt(tail(clim$Date,n=1))$year - as.POSIXlt(clim$Date[1])$year +1
lambda <- l/nyears   #average number of storms per year



################################
# generate under GAM modelling #
################################
library(gamlss)
faux<- function(x){ if(x<0){ y <-x+1}else{ y<-x}; y}
gaps <- sapply(ss.frac[2:l]-se.frac[1:(l-1)],faux)
tt <- ss.frac[1:(l-1)]
gam.fit <- gamlss(gaps ~ sin(2*pi*tt) + cos(2*pi*tt) + sin(4*pi*tt) + cos(4*pi*tt)
#                          + sin(8*pi*tt) + cos(8*pi*tt)
     , sigma.formula = ~  sin(2*pi*tt) + cos(2*pi*tt) + sin(4*pi*tt) + cos(4*pi*tt)
#                          + sin(8*pi*tt) + cos(8*pi*tt)
                , family="BE")

rGap <- function(){
  n.para   <- predictAll(gam.fit, newdata=data.frame(tt=yf), 
                           data=data.frame(gaps,tt))
  n.gap.aux<- qBE(runif(1), n.para$mu, n.para$sigma) #gap in year units
  as.integer(n.gap.aux*365*24) #gap in hours
}

rSind <- function(){ 
   sample(1:l,size=1)
}

Nstorms <- integer()
for(i in 1:n.s){ 
  gam.storms <- GenerateStormSample(rGap, rSind, pots.storm, clim, f.tide)
  Nstorms[i] <- nstorms
  write.csv(gam.storms, file=paste("Data.New/Storms-Gam-",i,".csv",sep=""), 
                                   row.names=F, quote= F)  
}
# test that the number of strorms per year is approximately the same
print(paste("Average number of stroms per year (simulation):", 
            mean(Nstorms)/50))
print(paste("Average number of stroms per year (original):", 
            lambda))

lambda <- mean(Nstorms)/50
# for the different clustering scenarios later on

##################################################################
# Repeat gam modelling introducing climate change sea level rise #
##################################################################

#read climate change
rcp26  <- read.table("./Data.Raw/rcp26_summid.txt")
rcp85  <- read.table("./Data.Raw/rcp85_summid.txt")

old.tide   <- f.tide$Level
iaux <- as.POSIXlt(f.tide$Date)$year-106
proj26 <- sapply(iaux,function(i){ rcp26$V2[i]})
proj85 <- sapply(iaux,function(i){ rcp85$V2[i]})
   
#scenario rcp26
f.tide$Level <- old.tide + proj26
for(i in 1:n.s){
  gam.storms <- GenerateStormSample(rGap, rSind, pots.storm, clim, f.tide)
  Nstorms[i] <- nstorms
  write.csv(gam.storms, file=paste("Data.New/Storms-rcp26-",i,".csv",sep=""),
                                   row.names=F, quote= F)
}
# test that the number of strorms per year is approximately the same
print(paste("Average number of stroms per year (simulation):",
            mean(Nstorms)/50))
print(paste("Average number of stroms per year (original):",
            lambda))



#scenario rcp85
f.tide$Level <- old.tide + proj85
gam.storms <- GenerateStormSample(rGap, rSind, pots.storm, clim, f.tide)

for(i in 1:n.s){
  gam.storms <- GenerateStormSample(rGap, rSind, pots.storm, clim, f.tide)
  Nstorms[i] <- nstorms
  write.csv(gam.storms, file=paste("Data.New/Storms-rcp85-",i,".csv",sep=""),
                                   row.names=F, quote= F)
}

# test that the number of strorms per year is approximately the same
print(paste("Average number of stroms per year (simulation):",
            mean(Nstorms)/50))
print(paste("Average number of stroms per year (original):",
            lambda))

#####################
# Equispaced storms in time
#####################
f.tide$Level <- old.tide
rGap   <- function(){ as.integer(365.25*24*1/lambda) }
rSind  <- function(){ sample(1:l,size=1) }
for(i in 1:n.s){
  eq.sp.storms <- GenerateStormSample(rGap, rSind, pots.storm, clim, f.tide)
  write.csv(eq.sp.storms, file=paste("Data.New/Storms-Equi-",i,".csv",sep=""),
                                     row.names=F, quote= F)
}

################################
# Poisson Distribution with the same number of storms per year
################################

rGap   <- function(){ as.integer(365.25*24*rexp(n=1,rate=lambda))}
rSind  <- function(){ sample(1:l,size=1) }
Nstorms <- integer()

for(i in 1:n.s){
  hom.pois.storms <- GenerateStormSample(rGap, rSind, pots.storm, clim, f.tide)
  Nstorms[i] <- nstorms
  write.csv(hom.pois.storms, file=paste("Data.New/Storms-Pois-",i,".csv",sep=""),
                                   row.names=F, quote= F)
}
# test that the number of strorms per year is approximately the same
print(paste("Average number of stroms per year (simulation):",
            mean(Nstorms)/50))
print(paste("Average number of stroms per year (original):",
            lambda))
