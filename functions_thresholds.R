combineSamples <- function(path, files, verbose = FALSE){
  combined.data <- NULL
  if(verbose == TRUE){
    cat("Reading Amplitude files.\n")
  }
  for(i in 1:length(files))
  { 
    nrlines <- NULL
    file.name <- file.path(path, files[i])
    if(length(count.fields(file.name)) > 0){
      sample.data <- read.table(file=file.name, header = TRUE, sep = ",")
      combined.data <- rbind(combined.data, sample.data )
    }
  }
  if(verbose == TRUE){
    cat("Amplitude files are combined.\n")
  }
  return(combined.data)
}
.makeBreaks <- function(min, max, breaks = 10){
  # not working with hist
  amplitude <- max - min
  amplitudeBreaks <- amplitude / breaks
  level <- min - 1
  result <- NULL
  for(i in 1:(breaks)){
    result <- c(result, level)
    level <- level + amplitudeBreaks
  }
  result <- c(result, (max+1))
  return(result)
}
maxAmplitude <- function(data, tData = NULL){
  result <- c(max(data[,1]), max(data[,2]))
  result <- thresholdData(tData = tData, amplitude = result, type = 'maxAmplitude')
  return(result)
}
minAmplitude <- function(data, tData = NULL){
  result <- c(min(data[,1]), min(data[,2]))
  result <- thresholdData(tData = tData, amplitude = result, type = 'minAmplitude')
  return(result)
}
meanCluster <- function(data, cluster = 1, channel = 1){
  result <- mean(data[$Cluster %in% cluster,channel])
  if(as.character(result) == "NaN"){result <- 0}
  return(result)
}
meanSdCluster <- function(data, cluster = 1, stdev = 1){
  if(cluster != 9){
    data <- data[data[,3] == cluster, ]
  }
  results <- c(mean(data[,1]), mean(data[,2]))
  results <- rbind(results, c((sd(data[,1])*stdev), (sd(data[,2]) * stdev)))
  colnames(results) <- c("Ch1","Ch2")
  rownames(results) <- c("mean", paste(stdev, "_sd", sep=""))
  return(results)
}
removeOutliers <- function(data, percentage = 1){
  percentage <- round((nrow(data)/100) * (percentage/2))
  selection <- c(1:percentage, (nrow(data)-percentage):nrow(data))
  selection <- !(1:nrow(data)  %in% selection)
  data <- data[order(data[,1]),]
  data <- data[selection,]
  selection <- c(1:percentage, (nrow(data)-percentage):nrow(data))
  selection <- !(1:nrow(data)  %in% selection)
  data <- data[order(data[,2]),]
  data <- data[selection,]
  return(data)
}
.thresholdHist <- function(x){
  x <- as.numeric(x)
  hist.data <- hist(x, breaks=15, plot=FALSE)
  hist.data <- rbind(hist.data$mids, hist.data$counts)
  hist.data <- hist.data[ ,-c(1:2,(dim(hist.data)[2] - 2):dim(hist.data)[2])]
  result <- mean(hist.data[1, hist.data[2,] == min(hist.data[2,])])
  return(result)
}
.thresholdKmeans2 <- function(x, nClusters = 2){
  x <- as.numeric(x)
  result <- NULL
  threshold <- kmeans(x = x, centers = nClusters)$centers
  if(dim(threshold)[1] == 2){result <- mean(threshold)}
  return(result)
}
.thresholdRanges <- function(x){
  x <- as.numeric(x)
  result <- NULL
  result <- (max(x) - min(x)) / 2
  return(result)
}
.thresholdKmeans4 <- function(x){
  # function is not optimized yet
  results <- kmeans(x[,1:2], 4)
  
  return(results$centers)
}
.thresholdDensity <- function(path = NULL, design = NULL, tData = NULL, breaks = 100, verbose = TRUE){
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
.thresholdDensityHist <- function(x, breaks = 100){
  hist.data <- hist(x, breaks = breaks, plot = FALSE)
  x <- hist.data$counts
  result <- NULL
  for(i in 1:(length(x)-5)){
    if(x[i+1] > x[i+2] & x[i+2] > x[i+3] &
       x[i+3] <= x[i+4]
    ) {
      result <- c(result, i+3)
    }
  }
  result <- result[1]
  results <- hist.data$mids[result]
  return(results)
}
setThresholdsManual <- function(ch1 = NULL, ch2 = NULL, tData = NULL, verbose = FALSE){
  if(is.null(ch1) == TRUE | is.null(ch2) == TRUE){
    stop("Thresholds not added correctly.\n")
  } else if (is.numeric(ch1) == FALSE | is.numeric(ch2) == FALSE){
    stop("Thresholds should be numeric.\n")
  } else if(verbose == TRUE){
    cat("Adding threshold manually.\n")
  }
  result <- c(ch1, ch2)
  result <- thresholdData(tData = tData, amplitude = result, type = 'threshold')
  return(result)
}

setThresholds <- function(data, algorithm = "densityhist", rm.outliers = TRUE, tData = NULL, verbose = FALSE){ 
  if(rm.outliers == TRUE){
    if("minOutlier" %in% row.names(tData) == TRUE) 
      { # min outliers
      data <- data[data[,1] > tData[row.names(tData) %in% "minOutlier", 1], ]
      data <- data[data[,2] > tData[row.names(tData) %in% "minOutlier", 2], ]
      if(verbose == TRUE){
        cat("\nMin outliers are removed from the data before setting thresholds.\n")
      }
    }
    if("maxOutlier" %in% row.names(tData) == TRUE) 
    { # max outliers
      data <- data[data[,1] < tData[row.names(tData) %in% "maxOutlier", 1], ]
      data <- data[data[,2] < tData[row.names(tData) %in% "maxOutlier", 2], ]
      if(verbose == TRUE){
        cat("Max outliers are removed from the data before setting thresholds.\n")
      }
    }
  }
  if(tolower(algorithm) == "hist" | tolower(algorithm) == "histogram"){
    if(verbose == TRUE){
      cat("\nSetting threshold based on 'histogram'.\n")
    }
    result <- c(.thresholdHist(x = data[,1]), .thresholdHist(x = data[,2]))
  }
  if(tolower(algorithm) == "ranges"){
    if(verbose == TRUE){
      cat("\nSetting threshold based on 'ranges'.\n")
    }
    result <- c(.thresholdRanges(x = data[,1]), .thresholdRanges(x = data[,2]))
  }
  if(tolower(algorithm) == "kmeans2"){
    if(verbose == TRUE){
      cat("\nSetting threshold based on 'kmeans' with 2 clusters per channel.\n")
    }
    result <- c(.thresholdKmeans2(x = data[,1]), .thresholdKmeans2(x = data[,2]))
  }
  if(tolower(algorithm) == "kmeans4"){
    if(verbose == TRUE){
      cat("\nSetting threshold based on 'kmeans' with 4 clusters for both channels.\n '1/3' of amplitude of Highest cluster will be set as threshold.\n")
    }
    result <- .thresholdKmeans4(x = data)$centers
    result <- c(max(result[,1])/3, max(result[,2]/3))
  }
  if(tolower(algorithm) == "densityhist"){
    if(verbose == TRUE){
      cat("\nSetting threshold based on 'densityhist'.\n")
    }
    result <- c(.thresholdDensityHist(x = data[,1]), .thresholdDensityHist(x = data[,2]))
  }
  
  result <- thresholdData(tData = tData, amplitude = result, type = 'threshold')
  return(result)
}
refineThresholdStDev <- function(data, stdev = 3, tData = NULL, verbose = FALSE){
  if (class(tData) == "matrix") 
  {
    if("threshold" %in% row.names(tData) == TRUE)
    {
      threshold <- tData[row.names(tData) %in% "threshold",]
    }else{
      stop("Threshold data is not given.\n")
      }
    cluster1 <- colSums(meanSdCluster(data, cluster = 1, stdev = stdev))
    cluster2 <- colSums(meanSdCluster(data, cluster = 2, stdev = stdev))
    cluster4 <- colSums(meanSdCluster(data, cluster = 4, stdev = stdev))
    
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
    result <- thresholdData(tData = tData, amplitude = threshold, type = 'threshold')
    return(result)
    
  }else{stop("Threshold data is not given as a matrix.\n")}
}
.createThresholdMatrix <- function(){
  # not used anymore? Remove?
  result <- matrix(data = NA, nrow = 7, ncol = 2, 
                   dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minAmplitude","minRain","maxRain"), c("ch1","ch2")))
  return(result)
}
.addThresholdData <- function(x = NULL, add, type = ''){
  # not used anymore? Remove?
  if (class(result != "matrix"))
  {
    result <- matrix(data = NA, nrow = 7, ncol = 2, 
                     dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minAmplitude","minRain","maxRain"), c("ch1","ch2")))
    
  }
  dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minAmplitude","minRain","maxRain"), c("ch1","ch2"))
  
}
thresholdData <- function(tData = NULL, amplitude = NULL, type = NULL) {
  if(class(type) != "NULL")
  {
    thresholds <- c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minAmplitude","mode","minRain","maxRain")
    if(tolower(type) %in% tolower(thresholds) == TRUE)
    {
      type <- thresholds[tolower(thresholds) %in% tolower(type)]
    } else {stop("No correct threshold type is given.\n")}
  }
  if(length(amplitude) == 2)
  {
    if(class(tData) == "matrix")
    {
      if(type %in% rownames(tData) == TRUE)
      {
        tData[grep(pattern = type, x = rownames(tData)),] <- amplitude
      }else
      {
        tData <- rbind(tData, type = amplitude)
        rownames(tData)[nrow(tData)] <- type
      }
    }else
    {
      tData <- matrix(data = amplitude, nrow = 1, ncol = 2, dimnames = list(c(type), c("ch1.Amplitude","ch2.Amplitude")))
    }
  }else{stop("Threshold data must have a length of 2.\n")}
  
  return(tData)
}
# [ ] - create new function that will create matrix or add settings to it.
# some fixed settings
# able to add own settings as well
# 