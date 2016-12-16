refineThresholdStDev <- function(x, stdev = 3, tData = NULL)
{
  if (class(tData) == "matrix") 
  {
    if("threshold" %in% row.names(tData) == TRUE)
    {
    threshold <- tData[row.names(tData) %in% "threshold",]
    }else{stop("Threshold data is not given.\n")}

    cluster1 <- colSums(calculateMeanSdCluster(x, cluster = 1, stdev = stdev))
    cluster2 <- colSums(calculateMeanSdCluster(x, cluster = 2, stdev = stdev))
    cluster4 <- colSums(calculateMeanSdCluster(x, cluster = 4, stdev = stdev))
    
    refined.channel1 <- max(c(cluster1[1], cluster4[1]), na.rm = TRUE)
    refined.channel2 <- max(c(cluster1[2], cluster2[2]), na.rm = TRUE)
    
    if(refined.channel1 < threshold[1])
    {
      threshold[1] <- mean(c(refined.channel1, threshold[1]))
    } else {cat("threshold for channel 1 could not be defined further.\n")}
    if(refined.channel2 < threshold[2])
    {
      threshold[2] <- mean(c(refined.channel2, threshold[2]))
    } else {cat("threshold for channel 2 could not be defined further.\n")}
    
    result <- thresholdData(tData = tData, amplitude = threshold, type = 'thresholdMeanStDev')
    return(result)
    
  }else{stop("Threshold data is not given as a matrix.\n")}
}