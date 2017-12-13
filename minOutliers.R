.removePercentage <- function(x, percentage = .5){
  x <- x[!(x %in% NA)]
  percentage <- round((length(x)/100) * (percentage/2))
  selection <- c(1:percentage, (length(x)-percentage):length(x))
  selection <- !(1:length(x)  %in% selection)
  x <- x[order(x)]
  x <- x[selection]
  return(x)
}
.removeOutliers <- function(x, percentage = .5){
  x <- x[!(x[,1] %in% NA),]
  percentage <- round((nrow(x)/100) * (percentage/2))
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
minOutliers <- function(data = NULL, verbose = TRUE){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  }
  if(verbose == TRUE){
    cat("Processing samples for outliers in the lower ranges.\n")
  }
  samples <- data@phenoData$sampleData['name',]
  for(i in 1:length(samples)){
    ### ANALYSIS CHANNEL 1
    channel1 <- as.numeric(data@assayData$Ch1.Amplitude[,i])
    hist.data <- hist(channel1, breaks=300, plot=FALSE)
    firstBigClusterCh1 <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
    
    if(TRUE %in% (hist.data$counts[1:firstBigClusterCh1] %in% 0 == TRUE)){
      closestZero <- max((1:firstBigClusterCh1)[hist.data$counts[1:firstBigClusterCh1] == 0])
      data@phenoData$ch1['minOutlier',i] <- hist.data$mids[closestZero]
    } else  {data@phenoData$ch1['minOutlier',i] <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) 
    if(data@phenoData$ch1['minOutlier',i] < 0){data@phenoData$ch1['minOutlier',i] <- 0 }
    }
    ### ANALYSIS CHANNEL 2
    channel2 <- as.numeric(data@assayData$Ch2.Amplitude[,i])
    hist.data <- hist(channel2, breaks=300, plot=FALSE)
    firstBigClusterCh2 <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
    
    if(TRUE %in% (hist.data$counts[1:firstBigClusterCh2] %in% 0 == TRUE)){
      closestZero <- max((1:firstBigClusterCh2)[hist.data$counts[1:firstBigClusterCh2] == 0])
      data@phenoData$ch2['minOutlier',i] <- hist.data$mids[closestZero]
    } else {data@phenoData$ch2['minOutlier',i] <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) 
    if(data@phenoData$ch2['minOutlier',i] < 0){data@phenoData$ch2['minOutlier',i] <- 0 }
    }
  }
  data <- .updateClusters(data)
  return(data)
}
