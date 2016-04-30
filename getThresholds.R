getThresholds <- function(x, algorithm = "hist", rm.outliers = TRUE)
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
  
  result <- matrix(data = result, nrow = 1, ncol = 2, dimnames = list(c("threshold"), c("ch1","ch2")))
  return(result)
}