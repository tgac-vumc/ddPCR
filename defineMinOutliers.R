defineMinOutliers <- function(x)
{
  channel1 <- as.numeric(x[,1])
  hist.data <- hist(channel1, breaks=200, plot=FALSE)
  firstBigCluster <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  closestZero <- max((1:firstBigCluster)[hist.data$counts[1:firstBigCluster] == 0])
  x[x[,1] < hist.data$mids[closestZero],3] <- 0
  
  channel2 <- as.numeric(x[,2])
  hist.data <- hist(channel2, breaks=200, plot=FALSE)
  firstBigCluster <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  closestZero <- max((1:firstBigCluster)[hist.data$counts[1:firstBigCluster] == 0])
  x[x[,2] < hist.data$mids[closestZero],3] <- 0
  
  return(x)
}

