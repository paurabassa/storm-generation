st      <- read.csv("New.Data/storms-Orig.csv")
st$Date <- as.POSIXct(strptime(st.org$Date, "%Y-%m-%d %H:%M:%S","GMT"))

st.dt <- abs(as.integer(difftime(st.org$Date[1:(l.org-1)],st.org$Date[2:l.org], 
                              units="hours")))

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
  while(st.dt[iaux]==1 & iaux <l.st){ clst <- clst+1; iaux <- iaux+1}
  if(iaux>=l.st) break
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

# convert start to date (add origin))
start <- as.POSIXct(start, origin="1970-01-01 00:00.00 UTC")
end   <- as.POSIXct(end,   origin="1970-01-01 00:00.00 UTC")
strm <- data.fame(yfst, yfend, start, end, dur, gaps, hs.m, res.m, fp.m, dir.av)

#proxi of the dependencies between variables
plot(data.frame(rank(strm$yfst,ties.method="r"),rank(strm$dur,ties.method="r"), rank(strm$gaps,ties.method="r"), rank(strm$hs.m,ties.method="r"), rank(strm$res.m,ties.method="r")))
