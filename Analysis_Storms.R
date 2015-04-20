# this file read a csv of storm occurence from New.Data and computes 
# some of the characteristics from Callaghan et. al

st      <- read.csv("Data.New/storms-Orig.csv")
#st      <- read.csv("New.Data/storms-Harm-2.csv")
st$Date <- as.POSIXct(strptime(st$Date, "%Y-%m-%d %H:%M:%S","GMT"))
l.org <-length(st$Date)

st.dt <- abs(as.integer(difftime(st$Date[1:(l.org-1)],st$Date[2:l.org], 
                              units="hours"))) #storm difference times"

source("R.functions/YearFraction.R")

iaux <- 1
l.st <- dim(st)[1]
n.st <- 1
dur   <- vector()
gaps  <- vector()
start <- vector()
end   <- vector()
yfst  <- vector()
yfend <- vector()
hs.m  <- vector()
fp.m  <- vector()
tm.m  <- vector()
res.m <- vector()
dir.av<- vector()

while(iaux < l.st){
  clst <-1
  #search for a jump in st Dates
  while(st.dt[iaux]==1 & iaux <l.st){ clst <- clst+1; iaux <- iaux+1} 
  if(iaux>=l.st) break
  # extract storm characteristics 
  dur[n.st]   <- clst
  gaps[n.st]  <- st.dt[iaux]
  start[n.st] <- st$Date[iaux-clst+1]
  end[n.st]   <- st$Date[iaux]
  yfst[n.st]  <- yearfraction( st$Date[iaux-clst+1]) 
  yfend[n.st] <- yearfraction( st$Date[iaux])
  hs.m[n.st]  <-  max(st$hs[(iaux-clst+1):iaux])
  res.m[n.st] <-  max(st$res[(iaux-clst+1):iaux])
  fp.m[n.st]  <-  max(st$fp[(iaux-clst+1):iaux])
  dir.av[n.st]<- mean(st$dir[(iaux-clst+1):iaux])
  n.st <- n.st + 1
  iaux <- iaux + 1
  print(iaux)
}

# convert start to date (add origin) and create a dataframe
start <- as.POSIXct(start, origin="1970-01-01 00:00.00 UTC")
end   <- as.POSIXct(end,   origin="1970-01-01 00:00.00 UTC")
strm  <- data.frame(yfst, yfend, start, end, dur, gaps, hs.m, res.m, fp.m, dir.av)

#scater plot of some variables
plot(strm[c("yfst", "dur", "gaps", "hs.m", "res.m", "fp.m", "dir.av")])

#proxi of the dependencies between variables
plot(data.frame(rank(strm$yfst,ties.method="r"),rank(strm$dur,ties.method="r"), rank(strm$gaps,ties.method="r"), rank(strm$hs.m,ties.method="r"), rank(strm$res.m,ties.method="r")))
