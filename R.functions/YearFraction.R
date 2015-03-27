yearfraction <- function(date){
# Function to compute the fraction of the year form a POSIX date 
# (takes into account leap years)
   aux   <- as.POSIXlt(date)
   year  <- aux$year
   lyear <- difftime(as.POSIXct(strptime(paste(1901+year,"-01-01",sep=""),
                                         "%Y-%m-%d")),
                     as.POSIXct(strptime(paste(1900+year,"-01-01",sep=""),
                                         "%Y-%m-%d")),
                     units="hours")
   fyear <- (24*aux$yday + aux$hour)/lyear[[1]]
}
