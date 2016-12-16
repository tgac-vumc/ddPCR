getMeanCluster <- function(x, cluster = 1, channel = 1)
{
  result <- mean(x[x$Cluster %in% cluster,channel])
  if(as.character(result) == "NaN"){result <- 0}
  return(result)
}