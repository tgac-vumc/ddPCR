getThresholdHist <- function(x)
{
  x <- as.numeric(x)
  hist.data <- hist(x, breaks=15, plot=FALSE)
  hist.data <- rbind(hist.data$mids, hist.data$counts)
  hist.data <- hist.data[ ,-c(1:2,(dim(hist.data)[2] - 2):dim(hist.data)[2])]
  result <- mean(hist.data[1, hist.data[2,] == min(hist.data[2,])])
  return(result)
}