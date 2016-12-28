getThresholdHist <- function(x){
  x <- as.numeric(x)
  hist.data <- hist(x, breaks=15, plot=FALSE)
  hist.data <- rbind(hist.data$mids, hist.data$counts)
  hist.data <- hist.data[ ,-c(1:2,(dim(hist.data)[2] - 2):dim(hist.data)[2])]
  result <- mean(hist.data[1, hist.data[2,] == min(hist.data[2,])])
  return(result)
}
getThresholdKmeans <- function(x, nClusters = 2){
  x <- as.numeric(x)
  result <- NULL
  threshold <- kmeans(x = x, centers = nClusters)$centers
  if(dim(threshold)[1] == 2){result <- mean(threshold)}
  return(result)
}
getThresholdRanges <- function(x){
  x <- as.numeric(x)
  result <- NULL
  result <- (max(x) - min(x)) / 2
  return(result)
}
getThresholdsKmeans <- function(x, rm.outliers = TRUE){
  if(rm.outliers == TRUE)
  {
    x <- x[!x[,3] == 0,]
  }
  results <- kmeans(x[,1:2], 4)
  return(results)
}
getThresholds <- function(x, algorithm = "hist", rm.outliers = TRUE, tData = NULL){ 
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
