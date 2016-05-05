defineClusters <- function(x, tData = NULL)
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
    
    if("minOutlier" %in% row.names(tData) == TRUE)
    {
      minOutlier <- tData[row.names(tData) %in% "minOutlier",]
      results[x[,1] < minOutlier[1] & x[,2] < minOutlier[2]] <- 0 # minOutlier : cluster 0
    }
    if("maxOutlier" %in% row.names(tData) == TRUE)
    {
      maxOutlier <- tData[row.names(tData) %in% "maxOutlier",]
      results[x[,1] > maxOutlier[1] & x[,2] > maxOutlier[2]] <- 0 # maxOutlier : cluster 0
    }
    
    x[,3] <- results
    }else{stop("Threshold data is not given as a matrix.\n")}
  return(x)
}