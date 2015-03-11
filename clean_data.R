#R file to clean the data and create a syncronized unique file with all 
#the relevant data. 

#read hindcast data
hcast <- read.csv("~/Post-Doc-LDN/Storm-Generation/Raw.Data/data.csv")
str(hcast)
hcast$Date <- as.POSIXct(strptime(hcast$Date, "%Y-%m-%d %H:%M:%S"))

#read past tides
res <- read.table("~/Post-Doc-LDN/Storm-Generation/Raw.Data/resid.dat")
M<- data.matrix(res)

tl <- length(res$V1)
v.res <- vector("numeric",tl*12)
for (i in 1:tl){ 
  for(j in 1:12){  v.res[12*(i-1)+j]=M[i,j] }
}
str(v.res)
plot(v.res, typ='l')

tot <- read.table("~/Post-Doc-LDN/Storm-Generation/Raw.Data/total.dat")
M<- data.matrix(tot)

tl <- length(tot$V1)
v.tot <- vector("numeric",tl*12)
for (i in 1:tl){ 
  for(j in 1:12){  v.tot[12*(i-1)+j]=M[i,j] }
}
str(v.tot)
plot(v.tot, typ='l')


# NOW put this into a data.frame format with the right date and time

tm <- as.POSIXct(strptime("0101199200", "%d%m%Y%H")) + 3600*(0:(length(v.res)-1))
tides <- data.frame(tm,v.tot,v.res, v.tot-v.res)
names(tides) <-c("Date", "t.tide", "res", "a.tide")


all <- merge(hcast, tides, by="Date")
o.clim <- all[,c("Date","hs", "fp", "dir", "U10", "V10", "t.tide", "res")]

plot(o.clim[o.clim$hs>2.5,c("hs", "fp", "dir", "res")])

storms.times <- as.POSIXlt(o.clim$Date[o.clim$hs>2.5])
plot(storms.times$yday,storms.times$year+1900,ylim=c(2013,1990))

rng <- range(storms.times$yday)
inc <- (rng[2]-rng[1])/12
seq(rng[1],rng[2],inc)
hist(storms.times$yday,breaks= seq(rng[1],rng[2],inc))


#read future tides

con <- file("~/Post-Doc-LDN/Storm-Generation/Raw.Data/6000591_2015-2064.txt")
lines <- readLines(con)
close(con)
which(nchar(lines)!=69)

lines2 <- lines[c(6:350645, 350651:438322)]

ff <- tempfile()
cat(file = ff, lines2, sep = "\n")
ftides <- read.fwf(ff, widths = c(18,10,10,10,12,12))    #
ftides$V1 <- as.POSIXct(strptime(as.character(ftides$V1), "%d/%m/%Y %H:%M"))
names(ftides)[1] <- "Date"
names(ftides)[2] <- "Level"
names(ftides)[3] <- "Speed"
names(ftides)[4] <- "Direc"
names(ftides)[5] <- "U-Comp"
names(ftides)[6] <- "V-Comp"


#library(copula)

