#
# R file to collect, clean and save relevant wave and tide data for its analysis
#
# several files are generated: 
# 
# 1 hiscast with the variables of interest every 3 hours f
#
#
 

#read wave hindcast data
hcast.all <- read.csv("~/Post-Doc-LDN/Storm-Generation/Raw.Data/data.csv")
#str(hcast)
hcast.all$Date <- as.POSIXct(strptime(hcast.all$Date, "%Y-%m-%d %H:%M:%S"))

waves <- hcast.all[!is.na(hcast.all$hs),c("Date","hs", "fp", "dir", "U10", "V10")]
#str(waves)

# this data is every 3 hours at the begining and the every hour at some point
# the following lines is to split the data into two files with uniform timescale
t1 <- min(which(as.logical(as.POSIXlt(waves[,1])$hour %%3)))
t.hourly <-(t1-1):length(waves$Date)
laux <- length(t.hourly)/3
t.3hours <- c(1:(t1-1),(t1-1)+3*(1:(laux-1)))
#which(as.logical(as.POSIXlt(hcast[t.3hours,1])$hour %%3)) 


write.csv(waves, file="./Clean.Data/waves-all.csv", row.names=F, quote=F)
write.csv(waves[t.hourly,], file="./Clean.Data/waves-hourly.csv", row.names=F, quote=F)
write.csv(waves[t.3hours,], file="./Clean.Data/waves-3hours.csv", row.names=F, quote=F)


#
#read residual tides and put the values into a vector. 
res   <- read.table("~/Post-Doc-LDN/Storm-Generation/Raw.Data/resid.dat")
M     <- data.matrix(res)
tl    <- length(res$V1)
v.res <- vector("numeric",tl*12)
for (i in 1:tl){ 
  for(j in 1:12){  v.res[12*(i-1)+j]=M[i,j] }
}
#str(v.res)
#plot(v.res, typ='l')

#read total tides and put the values into a vector. 
tot   <- read.table("~/Post-Doc-LDN/Storm-Generation/Raw.Data/total.dat")
M     <- data.matrix(tot)
tl    <- length(tot$V1)
v.tot <- vector("numeric",tl*12)
for (i in 1:tl){ 
  for(j in 1:12){  v.tot[12*(i-1)+j]=M[i,j] }
}
#str(v.tot)
#plot(v.tot, typ='l')

# Create a data.frame of tides with the corresponding time 
tm <- as.POSIXct(strptime("0101199200", "%d%m%Y%H")) + 3600*(0:(length(v.res)-1))
tides <- data.frame(tm, v.tot, v.res, v.tot-v.res)
names(tides) <-c("Date", "t.tide", "res", "a.tide")
write.csv(tides, file="./Clean.Data/tides-all.csv", row.names=F, quote=F)

#
# combine both files
#
climate <- merge(waves, tides, by="Date")

# the resulting merged set is also sampled every 3 hours at the begining 
# and then every hour at some point
# the following lines is to split the data into two files with uniform timescale
t1 <- min(which(as.logical(as.POSIXlt(climate[,1])$hour %%3)))
t.hourly <-(t1-1):length(climate$Date)
laux <- length(t.hourly)/3
t.3hours <- c(1:(t1-1),(t1-1)+3*(1:(laux-1)))

write.csv(climate, file = "./Clean.Data/clim-all.csv", row.names=F, quote=F)
write.csv(climate[t.hourly,], file = "./Clean.Data/clim-hourly.csv", row.names=F, quote=F)
write.csv(climate[t.3hours,], file = "./Clean.Data/clim-3hours.csv", row.names=F, quote=F)

#
# 2nd part: Read future files and write it in a second file
#
# read future tides

con <- file("./Raw.Data/6000591_2015-2064.txt")
lines <- readLines(con)
close(con)
#which(nchar(lines)!=69)

lines2 <- lines[c(6:350645, 350651:438322)] #this is just because the original 
#file has a header that hasn't been removed when concatenating files

ff <- tempfile()
cat(file = ff, lines2, sep = "\n")
ftides <- read.fwf(ff, widths = c(18,10,10,10,12,12))    #
#str(ftides)
ftides$V1 <- as.POSIXct(strptime(as.character(ftides$V1), "%d/%m/%Y %H:%M"))
names(ftides)[1] <- "Date"
names(ftides)[2] <- "Level"
names(ftides)[3] <- "Speed"
names(ftides)[4] <- "Direc"
names(ftides)[5] <- "U-Comp"
names(ftides)[6] <- "V-Comp"

write.csv(ftides, file = "./Clean.Data/future-tides.csv",row.names=F, quote=F)
#library(copula)

