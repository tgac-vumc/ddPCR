refineThresholdStDev <- function(x, stdev = 3, tData = NULL)
{
  if (class(tData) == "matrix") 
  {
    if("thresholdMeanStDev" %in% row.names(tData) == TRUE)
    {
      threshold <- tData[row.names(tData) %in% "thresholdMeanStDev",]
    }else{stop("thresholdMeanStDev data is not given.\n")}
    
   
    #cluster1 <- colSums(calculateMeanSdCluster(x, cluster = 1, stdev = stdev))
    #cluster2 <- colSums(calculateMeanSdCluster(x, cluster = 2, stdev = stdev))
    #cluster4 <- colSums(calculateMeanSdCluster(x, cluster = 4, stdev = stdev))
    
    #refined.channel1 <- max(c(cluster1[1], cluster4[1]), na.rm = TRUE)
    #refined.channel2 <- max(c(cluster1[2], cluster2[2]), na.rm = TRUE)
    

    
    result <- thresholdData(tData = tData, amplitude = threshold, type = 'minRain')
    result <- thresholdData(tData = tData, amplitude = threshold, type = 'maxRain')
    return(result)
    
  }else{stop(".\n")}
}