calculateMeanSdCluster <- function(x, cluster = 1, stdev = 1)
{
  if(cluster != 9){
    x <- filter(x, Cluster == cluster)
  }
  results <- c(mean(x[,1]), mean(x[,2]))
  results <- rbind(results, c((sd(x[,1])*stdev), (sd(x[,2]) * stdev)))
  colnames(results) <- c("Ch1","Ch2")
  rownames(results) <- c("mean", paste(stdev, "_sd", sep=""))
  return(results)
}