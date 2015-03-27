storms1      <- read.csv("New.Data/storms-Harm.csv")
storms1$Date <- as.POSIXct(strptime(storms1$Date, "%Y-%m-%d %H:%M:%S","GMT"))
storms2      <- read.csv("New.Data/storms-Pois.csv")
storms2$Date <- as.POSIXct(strptime(storms2$Date, "%Y-%m-%d %H:%M:%S","GMT"))


op <- par(mfrow=c(1,2))
plot(storms1[c("Date","hs")], ylim=c(2.4,9))
plot(storms2[c("Date","hs")], ylim=c(2.4,9))

par(op)
