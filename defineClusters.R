defineClusters <- function(x, tData = NULL, outlier = TRUE, rain = TRUE)
{
  if (class(tData) == "matrix") 
    {
    if("threshold" %in% row.names(tData) == TRUE)
      {
      threshold <- tData[row.names(tData) %in% "threshold",]
      }
    if("thresholdMeanStDev" %in% row.names(tData) == TRUE)
      {
      threshold <- tData[row.names(tData) %in% "thresholdMeanStDev",]
    }
    
    results <- rep(NA, dim(x)[1])
    results[x[,1] < threshold[1] & x[,2] < threshold[2]] <- 1 # ch1-ch2- : cluster 1
    results[x[,1] > threshold[1] & x[,2] < threshold[2]] <- 2 # ch1+ch2- : cluster 2
    results[x[,1] > threshold[1] & x[,2] > threshold[2]] <- 3 # ch1+ch2+ : cluster 3
    results[x[,1] < threshold[1] & x[,2] > threshold[2]] <- 4 # ch1-ch2+ : cluster 4
    if (outlier == TRUE){
      if("minOutlier" %in% row.names(tData) == TRUE)
      {
        minOutlier <- tData[row.names(tData) %in% "minOutlier",]
        results[x[,1] < minOutlier[1] & x[,2] < threshold[2]] <- 0 
        results[x[,2] < minOutlier[2] & x[,1] < threshold[1]] <- 0 # minOutlier : cluster 0
      }
      if("maxOutlier" %in% row.names(tData) == TRUE)
      {
        maxOutlier <- tData[row.names(tData) %in% "maxOutlier",]
        results[x[,1] > maxOutlier[1]] <- 0
        results[x[,2] > maxOutlier[2]] <- 0 # maxOutlier : cluster 0
      }
    }
    if (rain == TRUE){
      if("minRain" %in% row.names(tData) == TRUE & "maxRain" %in% row.names(tData) == TRUE)
      {
        minRain <- tData[row.names(tData) %in% "minRain",]
        maxRain <- tData[row.names(tData) %in% "maxRain",]
        results[x[,1] > minRain[1] & x[,1] < maxRain[1] & (results == 1 | results == 2) ] <- 5 # Rain ch1 : cluster 5
        results[x[,2] > minRain[2] & x[,2] < maxRain[2] & (results == 1 | results == 4) ] <- 5 # Rain ch2 : cluster 5
      }
    }

    x[,3] <- results
    }else{stop("Threshold data is not given as a matrix.\n")}
  return(x)
}