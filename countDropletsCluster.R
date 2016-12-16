countDropletsCluster <- function(x, cluster = 1)
{
  result <- sum(x$Cluster %in% cluster)
  return(result)
}