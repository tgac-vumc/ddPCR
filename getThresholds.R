removeOutliers <- function(x, cutoff = 5){
  percentage <- round((nrow(x)/100) * (cutoff/2))
  selection <- c(1:percentage, (nrow(x)-percentage):nrow(x))
  selection <- !(1:nrow(x)  %in% selection)
  x <- x[order(x[,1]),]
  x <- x[selection,]
  selection <- c(1:percentage, (nrow(x)-percentage):nrow(x))
  selection <- !(1:nrow(x)  %in% selection)
  x <- x[order(x[,2]),]
  x <- x[selection,]
  return(x)
}
getNegatives <- function(design = NULL, verbose = TRUE){
  if(is.null(design) == TRUE){
    stop("No design file is given.\n")
  } else {
    negatives <- design$Type == "neg"
    if(TRUE %in% negatives == TRUE){
      design <- design[negatives,]
    } else {
      stop("No negative control was found.\n")
    }
  } 
  return(design)
}
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
getThresholdDensity <- function(path = NULL, design = NULL, tData = NULL, breaks = 100, verbose = TRUE){
  if(is.null(design) == TRUE){
    stop("No design file has been given.\n")
  } else {
    if(verbose == TRUE){
      cat("processing negative samples.\n")
    }
    design <- getNegatives(design)
    x <- combineSamples(path = path, files = design$File)
    x <- removeOutliers(x = x, cutoff = 3)
    maxAmplitude <- getMaxAmplitude(x)
    densityCh1 <- hist(x[,1], breaks=breaks, plot=FALSE)
    densityCh1 <- densityCh1$mids[densityCh1$density == max(densityCh1$density)]
    densityCh2 <- hist(x[,2], breaks=breaks, plot=FALSE)
    densityCh2 <- densityCh2$mids[densityCh2$density == max(densityCh2$density)]
    
    for(stdevSetting in 1:10){
      meanPlusSd <- colSums(calculateMeanSdCluster(x, stdev = stdevSetting, cluster = 9))
      if((meanPlusSd[1] > maxAmplitude[1] & meanPlusSd[2] > maxAmplitude[2]) == TRUE){
        stdevSetting <- (stdevSetting + 10)
        meanPlusSd <- colSums(calculateMeanSdCluster(x, stdev = stdevSetting, cluster = 9))
        tData <- thresholdData(tData = tData, amplitude = meanPlusSd, type = 'threshold')
        break()
      }
    }
    if(verbose == TRUE){
      cat("Threshold from negative is set by 'mean +", stdevSetting,"sd'.\n")
    }
  }
  return(tData)
} 
getThresholdDensityHist <- function(x, breaks = 100){
  hist.data <- hist(x, breaks = breaks, plot = FALSE)
  x <- hist.data$counts
  result <- NULL
  for(i in 1:(length(x)-5)){
    if(x[i+1] > x[i+2] & x[i+2] > x[i+3] &
       x[i+3] < x[i+4]
    ) {
      result <- c(result, i+3)
    }
  }
  result <- result[1]
  results <- hist.data$mids[result]
  return(results)
}

getThresholds <- function(x, algorithm = "densityhist", rm.outliers = TRUE, tData = NULL, verbose = FALSE){ 
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
    if(verbose == TRUE){
      cat("Setting threshold based on 'histogram'.\n")
    }
    result <- c(getThresholdHist(x = x[,1]), getThresholdHist(x = x[,2]))
  }
  if(tolower(algorithm) == "ranges")
  {
    if(verbose == TRUE){
      cat("Setting threshold based on 'ranges'.\n")
    }
    result <- c(getThresholdRanges(x = x[,1]), getThresholdRanges(x = x[,2]))
  }
  if(tolower(algorithm) == "kmeans")
  {
    if(verbose == TRUE){
      cat("Setting threshold based on 'kmeans'.\n")
    }
    result <- c(getThresholdRanges(x = x[,1]), getThresholdRanges(x = x[,2]))
  }
  if(tolower(algorithm) == "densityhist")
  {
    if(verbose == TRUE){
      cat("Setting threshold based on 'densityhist'.\n")
    }
    result <- c(getThresholdDensityHist(x = x[,1]), getThresholdDensityHist(x = x[,2]))
  }
  
  result <- thresholdData(tData = tData, amplitude = result, type = 'threshold')
  return(result)
}
refineThresholdStDev <- function(x, stdev = 3, tData = NULL, verbose = FALSE){
  if (class(tData) == "matrix") 
  {
    if("threshold" %in% row.names(tData) == TRUE)
    {
      threshold <- tData[row.names(tData) %in% "threshold",]
    }else{
      stop("Threshold data is not given.\n")
      }
    cluster1 <- colSums(calculateMeanSdCluster(x, cluster = 1, stdev = stdev))
    cluster2 <- colSums(calculateMeanSdCluster(x, cluster = 2, stdev = stdev))
    cluster4 <- colSums(calculateMeanSdCluster(x, cluster = 4, stdev = stdev))
    
    refined.channel1 <- max(c(cluster1[1], cluster4[1]), na.rm = TRUE)
    refined.channel2 <- max(c(cluster1[2], cluster2[2]), na.rm = TRUE)
    
    if(refined.channel1 < threshold[1])
    {
      if(verbose == TRUE){
        cat("Refining threshold with Mean and Standard Deviation setting of",stdev, " for channel 1.\n")
      }
      threshold[1] <- mean(c(refined.channel1, threshold[1]))
    } else {
      if(verbose == TRUE){
      cat("threshold for channel 1 could not be defined further.\n")
      }
    }
    if(refined.channel2 < threshold[2])
    {
      if(verbose == TRUE){
        cat("Refining threshold with Mean and Standard Deviation setting of",stdev, " for channel 2.\n")
      }
      threshold[2] <- mean(c(refined.channel2, threshold[2]))
    } else {
      if(verbose == TRUE){
        cat("threshold for channel 2 could not be defined further.\n")
        }
    }
    result <- thresholdData(tData = tData, amplitude = threshold, type = 'thresholdMeanStDev')
    return(result)
    
  }else{stop("Threshold data is not given as a matrix.\n")}
}
