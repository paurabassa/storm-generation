# 
#  R script to compare the data generated with the two different models 
#  with the original set of data. 
#

#read data
st.org      <- read.csv("New.Data/storms-Orig.csv")
st.org$Date <- as.POSIXct(strptime(st.org$Date, "%Y-%m-%d %H:%M:%S","GMT"))
st.hrm      <- read.csv("New.Data/storms-Harm-2.csv")
st.hrm$Date <- as.POSIXct(strptime(st.hrm$Date, "%Y-%m-%d %H:%M:%S","GMT"))
st.gam      <- read.csv("New.Data/storms-Pois-2.csv")
st.gam$Date <- as.POSIXct(strptime(st.gam$Date, "%Y-%m-%d %H:%M:%S","GMT"))

# plot Date vs hs for the different data sets
op <- par(mfrow=c(1,3))
plot(st.org[c("Date","hs")], ylim=c(2.4,7))
plot(st.hrm[c("Date","hs")], ylim=c(2.4,7))
plot(st.gam[c("Date","hs")], ylim=c(2.4,7))
par(op)

#quantile-quantile plots of hs
op <- par(mfrow=c(1,2))
qqplot(st.org$hs, st.hrm$hs, xlim=c(2.5,6.5), ylim=c(2.5,6.5))
abline(0,1)
qqplot(st.org$hs, st.gam$hs, xlim=c(2.5,6.5), ylim=c(2.5,6.5))
abline(0,1)
par(op)

savePlot("QQ-hs.png",type="png")

#quantile-quantile plots of interarrival times 
l.org <- length(st.org$Date)
l.hrm <- length(st.hrm$Date)
l.gam <- length(st.gam$Date)

dt.org <- abs(as.integer(difftime(st.org$Date[1:(l.org-1)],st.org$Date[2:l.org], units="hours")))
dt.hrm <- abs(as.integer(difftime(st.hrm$Date[1:(l.hrm-1)],st.hrm$Date[2:l.hrm], units="hours")))
dt.gam <- abs(as.integer(difftime(st.gam$Date[1:(l.gam-1)],st.gam$Date[2:l.gam], units="hours")))

gaps.org <- dt.org[dt.org>1]
gaps.hrm <- dt.hrm[dt.hrm>1]
gaps.gam <- dt.gam[dt.gam>1]

op <- par(mfrow=c(1,2))
qqplot(gaps.org, gaps.hrm, xlim=c(0,8000),  ylim=c(0,8000))
abline(0,1)
qqplot(gaps.org, gaps.gam, xlim=c(0,8000),  ylim=c(0,8000))
abline(0,1)
par(op)

savePlot("QQ-interarr.png",type="png")

#
# Rest of variables are direct resample of the origial data in both cases, so 
# they should look pretty similar 
#
#quantile-quantile plots of hs
op <- par(mfrow=c(1,2))
lims<-range(st.org$fp, st.hrm$fp, st.gam$fp )
qqplot(st.org$fp, st.hrm$fp,xlim=lims, ylim=lims)
abline(0,1)
qqplot(st.org$fp, st.gam$fp,xlim=lims, ylim=lims)
abline(0,1)
par(op)

