getThresholds <- function(x, algorithm = "hist", rm.outliers = TRUE, tData = NULL)
{ 
  if(rm.outliers == TRUE)
  {
    x <- x[!x[,3] == 0,]
  }
  if(tolower(algorithm) == "hist" | tolower(algorithm) == "histogram")
  {
    result <- c(getThresholdHist(x = x[,1]), getThresholdHist(x = x[,2]))
  }
  if(tolower(algorithm) == "ranges")
  {
    result <- c(getThresholdRanges(x = x[,1]), getThresholdRanges(x = x[,2]))
  }
  if(tolower(algorithm) == "kmeans")
  {
    result <- c(getThresholdRanges(x = x[,1]), getThresholdRanges(x = x[,2]))
  }
  
  result <- thresholdData(tData = tData, amplitude = result, type = 'threshold')
  return(result)
}