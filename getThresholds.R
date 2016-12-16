getThresholds <- function(x, algorithm = "hist", rm.outliers = TRUE, tData = NULL)
{ 
  if(rm.outliers == TRUE)
  {
    if("minOutlier" %in% row.names(tData) == TRUE) 
      { # min outliers
      x <- x[x[,1] > tData[row.names(tData) %in% "minOutlier", 1], ]
      x <- x[x[,2] > tData[row.names(tData) %in% "minOutlier", 2], ]
    }
    if("maxOutlier" %in% row.names(tData) == TRUE) 
    { # max outliers
      x <- x[x[,1] < tData[row.names(tData) %in% "maxOutlier", 1], ]
      x <- x[x[,2] < tData[row.names(tData) %in% "maxOutlier", 2], ]
    }
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