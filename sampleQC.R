sampleQC <- function(x, sample.type="")
{ # input = amplitide data with defined clusters
  if(countDropletsCluster(x = x, cluster = c(1,2,3,4)) < 10000)
  {
    result <- "CHECK DROPLET COUNT"
  } else  {
    result <- "DROPLET COUNT OK"
  }
  if(tolower(sample.type) == "pos")
  {
    if(countDropletsCluster(x = x, cluster=1) == 0){result <- c(result,"NO DROPLETS CLUSTER 1")}
    if(countDropletsCluster(x = x, cluster=2) == 0){result <- c(result,"NO DROPLETS CLUSTER 2")}
    if(countDropletsCluster(x = x, cluster=3) == 0){result <- c(result,"NO DROPLETS CLUSTER 3")}
    if(countDropletsCluster(x = x, cluster=4) == 0){result <- c(result,"NO DROPLETS CLUSTER 4")}
  }
  if(tolower(sample.type) == "ntc")
  {
    if(countDropletsCluster(x = x, cluster = c(2,3,4)) > 0){result <- c(result,"FALSE POSITIVE FOUND")}
  }
  result <- paste(result, collapse=": ")
  return(result)
}