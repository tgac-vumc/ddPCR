defineMinOutliers <- function(x, tData = NULL)
{
  channel1 <- as.numeric(x[,1])
  hist.data <- hist(channel1, breaks=300, plot=FALSE)
  firstBigCluster <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  if(TRUE %in% (hist.data$counts[1:firstBigCluster] %in% 0 == TRUE))
    {
    closestZero <- max((1:firstBigCluster)[hist.data$counts[1:firstBigCluster] == 0])
    minOutlierCh1 <- hist.data$mids[closestZero]
    } else { minOutlierCh1 <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) }
  
  channel2 <- as.numeric(x[,2])
  hist.data <- hist(channel2, breaks=300, plot=FALSE)
  firstBigCluster <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  if(TRUE %in% (hist.data$counts[1:firstBigCluster] %in% 0 == TRUE))
  {
    closestZero <- max((1:firstBigCluster)[hist.data$counts[1:firstBigCluster] == 0])
    minOutlierCh2 <- hist.data$mids[closestZero]
  } else { minOutlierCh2 <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) }
  
  results <- thresholdData(tData = tData, amplitude = c(minOutlierCh1, minOutlierCh2), type = 'minOutlier')
  
  return(results)
}

