ExtractPOTPointProcess <- function(ts, threshold){
# Given a time series and a threshold this function computes the point 
# process of excedants. It returns a list with 3 components
# p.exd  ->  position excedance
# c.siz  ->  cluster size                (for storms this is duration)
# n.exc  ->  distance to next excedance  (for storms this is inter arr. time)

  output   <- list(p.exc = vector(mode="numeric"),
                   c.siz = vector(mode="numeric"),
                   nex.exc = vector(mode="numeric"))
  booleans <- ts > threshold
  i.time <- 1
  start.times<- vector()
  n.exc <- 0
  while(i.time < length(booleans)){
    while(!booleans[i.time] && i.time < length(booleans))
      i.time <- i.time+1   #find the first excedance 
    n.exc <- n.exc + 1
#    print("excedance number= ")
#    print(n.exc)

#    print("excedance at pos= ")
#    print(i.time)
    output$p.exc[n.exc] <- i.time

    clust.size <- 0
    while(booleans[i.time] && i.time < length(booleans)){
      i.time     <- i.time+1
      clust.size <- clust.size + 1
    }
#    print("cluster size =")
#    print(clust.size)
    output$c.siz[n.exc] <- clust.size
  }
  for(i in 1:(n.exc-1))
    output$nex.exc[i] <- output$p.exc[i+1] -
                               (output$p.exc[i] + output$c.siz[i])
  output
}

