defineMinOutliers <- function(x)
{
  channel1 <- as.numeric(x[,1])
  hist.data <- hist(channel1, breaks=300, plot=FALSE)
  firstBigCluster <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  closestZero <- max((1:firstBigCluster)[hist.data$counts[1:firstBigCluster] == 0])
  MinOutlierCh1 <- hist.data$mids[closestZero]
  
  channel2 <- as.numeric(x[,2])
  hist.data <- hist(channel2, breaks=300, plot=FALSE)
  firstBigCluster <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  closestZero <- max((1:firstBigCluster)[hist.data$counts[1:firstBigCluster] == 0])
  MinOutlierCh2 <- hist.data$mids[closestZero]
  
  results <- c(MinOutlierCh1 = MinOutlierCh1, MinOutlierCh2 = MinOutlierCh2)
  
  return(results)
}

