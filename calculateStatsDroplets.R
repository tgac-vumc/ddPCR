calculateStatsDroplets <- function(x)
{ # input = amplitide data with defined clusters
  col.names <- c("Positives", "Negatives", "Ch1-Ch2-", "Ch1+Ch2-", "Ch1+Ch2+", "Ch1-Ch2+", "Outlier", "AcceptedDroplets")
  results <- matrix(0, nrow=2, ncol=length(col.names), dimnames = list(NULL, col.names))
  results[1, colnames(results) == "Positives"] <- countDropletsCluster(x, c(2,3)) #Positives
  results[2, colnames(results) == "Positives"] <- countDropletsCluster(x, c(3,4)) #Positives
  results[1, colnames(results) == "Negatives"] <- countDropletsCluster(x, c(1,4))# Negatives
  results[2, colnames(results) == "Negatives"] <- countDropletsCluster(x, c(1,2)) # Negatives
  results[1:2, colnames(results) == "Ch1-Ch2-"] <- countDropletsCluster(x, 1)
  results[1:2, colnames(results) == "Ch1+Ch2-"] <- countDropletsCluster(x, 2)
  results[1:2, colnames(results) == "Ch1+Ch2+"] <- countDropletsCluster(x, 3)
  results[1:2, colnames(results) == "Ch1-Ch2+"] <- countDropletsCluster(x, 4)
  results[1:2, colnames(results) == "Outlier"] <- countDropletsCluster(x, 0)
  results[1:2, colnames(results) == "AcceptedDroplets"] <- countDropletsCluster(x, c(1,2,3,4)) #AcceptedDroplets
  return(results)
}