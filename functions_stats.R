calculateStatsAmplitudes <- function(x){
  col.names <- c("MeanAmplitudeofPositives", "MeanAmplitudeofNegatives", "MeanAmplitudeTotal")
  results <- matrix(0, nrow = 2, ncol = length(col.names), dimnames = list(NULL, col.names))
  results[1, colnames(results) == "MeanAmplitudeofPositives"] <- round(getMeanCluster(x, cluster=c(2,3),1), digits = 2)
  results[2, colnames(results) == "MeanAmplitudeofPositives"] <- round(getMeanCluster(x, cluster=c(4,3),2), digits = 2)
  results[1, colnames(results) == "MeanAmplitudeofNegatives"] <- round(getMeanCluster(x, cluster=c(1,4),1), digits = 2)
  results[2, colnames(results) == "MeanAmplitudeofNegatives"] <- round(getMeanCluster(x, cluster=c(1,2),2), digits = 2)
  results[1, colnames(results) == "MeanAmplitudeTotal"] <- round(getMeanCluster(x, cluster=c(1,2,3,4),1), digits = 2)
  results[2, colnames(results) == "MeanAmplitudeTotal"] <- round(getMeanCluster(x, cluster=c(1,2,3,4),2), digits = 2)
  return(results)
}
CalculateStatsCopies <- function(x, interations = 10){ # input = amplitide data with defined clusters
  col.names <- c("CopiesPer1ul","CopiesPer20ulWell")
  results <- matrix(0, nrow = 2, ncol = length(col.names), dimnames = list(NULL, col.names))
  channel.1 <- countDropletsCluster(x, c(2,3)) #Positives channel 1
  channel.2 <- countDropletsCluster(x, c(3,4)) #Positives channel 1
  total <- countDropletsCluster(x, c(1,2,3,4)) #AcceptedDroplets
  results[1,colnames(results) == "CopiesPer1ul"] <- round(calculateCopies(posCount = channel.1, count = total, vDroplet = 0.85), digits = 1)
  results[2,colnames(results) == "CopiesPer1ul"] <- round(calculateCopies(posCount = channel.2, count = total, vDroplet = 0.85), digits = 1)
  results[1:2,colnames(results) == "CopiesPer20ulWell"] <- results[1:2,colnames(results) == "CopiesPer1ul"]*20
  return(results)
}
calculateStatsDroplets <- function(x){ # input = amplitide data with defined clusters
  col.names <- c("Positives", "Negatives", "Ch1-Ch2-", "Ch1+Ch2-", "Ch1+Ch2+", "Ch1-Ch2+", "Outlier", "Rain", "AcceptedDroplets")
  results <- matrix(0, nrow=2, ncol=length(col.names), dimnames = list(NULL, col.names))
  results[1, colnames(results) == "Positives"]  <- countDropletsCluster(x, c(2,3)) #Positives
  results[2, colnames(results) == "Positives"]  <- countDropletsCluster(x, c(3,4)) #Positives
  results[1, colnames(results) == "Negatives"]  <- countDropletsCluster(x, c(1,4))# Negatives
  results[2, colnames(results) == "Negatives"]  <- countDropletsCluster(x, c(1,2)) # Negatives
  results[1:2, colnames(results) == "Ch1-Ch2-"] <- countDropletsCluster(x, 1)
  results[1:2, colnames(results) == "Ch1+Ch2-"] <- countDropletsCluster(x, 2)
  results[1:2, colnames(results) == "Ch1+Ch2+"] <- countDropletsCluster(x, 3)
  results[1:2, colnames(results) == "Ch1-Ch2+"] <- countDropletsCluster(x, 4)
  results[1:2, colnames(results) == "Outlier"]  <- countDropletsCluster(x, 0)
  results[1:2, colnames(results) == "Rain"]     <- countDropletsCluster(x, 5)
  results[1:2, colnames(results) == "AcceptedDroplets"] <- countDropletsCluster(x, c(1,2,3,4)) #AcceptedDroplets
  return(results)
}
calculateStatsRatioFract <- function(x){ # input = output get.statistics.copies()
  if ("CopiesPer1ul" %in% colnames(x) == TRUE)
  {
    col.names <- c("Ratio", "FractionalAbundance")
    results <- matrix(0, nrow = 2, ncol = length(col.names), dimnames = list(NULL, col.names))
    results[1:2, colnames(results) == "Ratio"] <- x[1, colnames(x) == "CopiesPer1ul"] / x[2, colnames(x) == "CopiesPer1ul"]
    results[2, colnames(results) == "Ratio"] <- x[2, colnames(x) == "CopiesPer1ul"] / x[1, colnames(x) == "CopiesPer1ul"]
    results[1, colnames(results) == "FractionalAbundance"] <- x[1, colnames(x) == "CopiesPer1ul"] / (x[1, colnames(x) == "CopiesPer1ul"]+x[2, colnames(x) == "CopiesPer1ul"]) * 100
    results[2, colnames(results) == "FractionalAbundance"] <- x[2, colnames(x) == "CopiesPer1ul"] / (x[1, colnames(x) == "CopiesPer1ul"]+x[2, colnames(x) == "CopiesPer1ul"]) * 100
    return(results)
  }
}
calculateCopies <-  function(posCount, count, vDroplet=0.91, volume=1){ # concentration in copies / user defined volume
  negCount <- count-posCount
  result <- ((-log(negCount/count)/vDroplet)) * 1000 * volume
  return(result)
}
calculateStatsDroplets <- function(x){ # input = amplitide data with defined clusters
  col.names <- c("Positives", "Negatives", "Ch1-Ch2-", "Ch1+Ch2-", "Ch1+Ch2+", "Ch1-Ch2+", "Outlier", "Rain", "AcceptedDroplets")
  results <- matrix(0, nrow=2, ncol=length(col.names), dimnames = list(NULL, col.names))
  results[1, colnames(results) == "Positives"]  <- countDropletsCluster(x, c(2,3)) #Positives
  results[2, colnames(results) == "Positives"]  <- countDropletsCluster(x, c(3,4)) #Positives
  results[1, colnames(results) == "Negatives"]  <- countDropletsCluster(x, c(1,4))# Negatives
  results[2, colnames(results) == "Negatives"]  <- countDropletsCluster(x, c(1,2)) # Negatives
  results[1:2, colnames(results) == "Ch1-Ch2-"] <- countDropletsCluster(x, 1)
  results[1:2, colnames(results) == "Ch1+Ch2-"] <- countDropletsCluster(x, 2)
  results[1:2, colnames(results) == "Ch1+Ch2+"] <- countDropletsCluster(x, 3)
  results[1:2, colnames(results) == "Ch1-Ch2+"] <- countDropletsCluster(x, 4)
  results[1:2, colnames(results) == "Outlier"]  <- countDropletsCluster(x, 0)
  results[1:2, colnames(results) == "Rain"]     <- countDropletsCluster(x, 5)
  results[1:2, colnames(results) == "AcceptedDroplets"] <- countDropletsCluster(x, c(1,2,3,4)) #AcceptedDroplets
  return(results)
}
convertCopiesToNanogram <- function(x){
  x <- as.numeric(x) / 15.152
  return(x)
}
countDropletsCluster <- function(data, cluster = 1){
  result <- sum(x$Cluster %in% cluster)
  return(result)
}
dropletCountText <- function(x){
  results <- NULL
  results <- list(clusters = c(cluster.1 = sum(x == 1), cluster.2 = sum(x == 2), cluster.3 = sum(x == 3), cluster.4 = sum(x == 4), outlier = sum(x == 0), rain = sum(x == 5)))
  results <- paste("Ch1-Ch2-:", results$clusters[1],
                   "   Ch1+Ch2-:", results$clusters[2],
                   "   Ch1+Ch2+:", results$clusters[3],
                   "   Ch1-Ch2+:", results$clusters[4], 
                   "   outlier:", results$clusters[5], 
                   "   rain:", results$clusters[6], 
                   sep="")
  return(results)
}