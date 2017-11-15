.defineClusters <- function(data, tData = NULL, outlier = TRUE, rain = TRUE){
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
    if (outlier == TRUE){
      if("minOutlier" %in% row.names(tData) == TRUE)
      {
        minOutlier <- tData[row.names(tData) %in% "minOutlier",]
        results[x[,1] < minOutlier[1] & x[,2] < threshold[2]] <- 0 
        results[x[,2] < minOutlier[2] & x[,1] < threshold[1]] <- 0 # minOutlier : cluster 0
      }
      if("maxOutlier" %in% row.names(tData) == TRUE)
      {
        maxOutlier <- tData[row.names(tData) %in% "maxOutlier",]
        results[x[,1] > maxOutlier[1]] <- 0
        results[x[,2] > maxOutlier[2]] <- 0 # maxOutlier : cluster 0
      }
    }
    if (rain == TRUE){
      if("minRain" %in% row.names(tData) == TRUE & "maxRain" %in% row.names(tData) == TRUE)
      {
        minRain <- tData[row.names(tData) %in% "minRain",]
        maxRain <- tData[row.names(tData) %in% "maxRain",]
        results[x[,1] > minRain[1] & x[,1] < maxRain[1] & (results == 1 | results == 2) ] <- 5 # Rain ch1 : cluster 5
        results[x[,2] > minRain[2] & x[,2] < maxRain[2] & (results == 1 | results == 4) ] <- 5 # Rain ch2 : cluster 5
      }
    }

    x[,3] <- results
    }else{stop("Threshold data is not given as a matrix.\n")}
  return(x)
}
.defineColor <- function(x, density = NULL){
  result <- rep(NA, length(x))
  result[grep(pattern = 0, x = x)] <- "#FF0000"
  result[grep(pattern = 1, x = x)] <- "#000000"
  result[grep(pattern = 2, x = x)] <- "#0033FF" 
  result[grep(pattern = 3, x = x)] <- "#FF6600"
  result[grep(pattern = 4, x = x)] <- "#00CC00"
  result[grep(pattern = 5, x = x)] <- "#FF00CC"
  if(class(density) != "NULL")
  {
    result <- paste(result, as.character(density), sep = "") 
  }
  return(result)
}
defineMinOutliers <- function(data, tData = NULL){
  channel1 <- as.numeric(x[,1])
  hist.data <- hist(channel1, breaks=300, plot=FALSE)
  firstBigClusterCh1 <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  if(TRUE %in% (hist.data$counts[1:firstBigClusterCh1] %in% 0 == TRUE))
  {
    closestZero <- max((1:firstBigClusterCh1)[hist.data$counts[1:firstBigClusterCh1] == 0])
    minOutlierCh1 <- hist.data$mids[closestZero]
  } else  { minOutlierCh1 <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) 
  if(minOutlierCh1 < 0){ minOutlierCh1 <- 0 }
  }
  
  channel2 <- as.numeric(x[,2])
  hist.data <- hist(channel2, breaks=300, plot=FALSE)
  firstBigClusterCh2 <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
  if(TRUE %in% (hist.data$counts[1:firstBigClusterCh2] %in% 0 == TRUE))
  {
    closestZero <- max((1:firstBigClusterCh2)[hist.data$counts[1:firstBigClusterCh2] == 0])
    minOutlierCh2 <- hist.data$mids[closestZero]
  } else { minOutlierCh2 <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) 
  if(minOutlierCh2 < 0){ minOutlierCh2 <- 0 }
  }
  
  results <- thresholdData(tData = tData, amplitude = c(minOutlierCh1, minOutlierCh2), type = 'minOutlier')
  
  return(results)
}
defineTheRain <- function(data, stdev = c(3,3), tData = NULL, verbose = TRUE){
  if (class(tData) == "matrix") {
    if("thresholdMeanStDev" %in% row.names(tData) == TRUE){
      threshold <- tData[row.names(tData) %in% "thresholdMeanStDev",]
    }else if (class(tData) == "matrix") {
      if("threshold" %in% row.names(tData) == TRUE){
        threshold <- tData[row.names(tData) %in% "threshold",]
      }else {stop("Threshold data is not given.\n")}
    }
    if(length(stdev) == 1){
      stdevLow <- stdev[1]
      stdevHigh <- stdev[1]
    }else if(length(stdev) == 2){
      stdevLow <- stdev[1]
      stdevHigh <- stdev[2]
    }else {stop("Standard deviation setting to define the Rain cut-off is not in the correct format.\n")}
    
    if(verbose == TRUE){
      cat("DefineTheRain is calculated with a Standard Deviation of",stdev[1], "for ch1 and", stdev[2], "for ch2.\n")
    }  
    
    x <- defineClusters(x = x, tData = tData)
    rainResult <- matrix(data = NA, nrow = 2, ncol = 2, 
                         dimnames = list(c("minRain", "maxRain"), c("ch1.Amplitude","ch2.Amplitude")), 
                         byrow = TRUE)
    ch1Low <- calculateMeanSdCluster(x, cluster = c(1), stdev = stdevLow)
    rainResult[1,1] <- ch1Low[1,1] + ch1Low[2,1]
    ch1High <- calculateMeanSdCluster(x, cluster = c(2), stdev = stdevHigh)
    rainResult[2,1] <- ch1High[1,1] - ch1High[2,1]
    ch2Low <- calculateMeanSdCluster(x, cluster = c(1), stdev = stdevLow)
    rainResult[1,2] <- ch2Low[1,2] + ch2Low[2,2]
    ch2High <- calculateMeanSdCluster(x, cluster = c(4), stdev = stdevHigh)
    rainResult[2,2] <- ch2High[1,2] - ch2High[2,2]
    
    tData <- thresholdData(tData = tData, amplitude = rainResult[1,], type = 'minRain')
    tData <- thresholdData(tData = tData, amplitude = rainResult[2,], type = 'maxRain')
    return(tData)
  }
}
