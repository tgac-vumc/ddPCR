.removeOutliers <- function(x,y, percentage = 1){
  if(length(x) == length(y)){
    x <- cbind(x, y)
  }
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
.thresholdHist <- function(x, breaks = 15, strict = FALSE){
  x <- x[!(x %in% NA)]
  if(strict == TRUE){
    breaks <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
  }
  hist.data <- hist(x, breaks = breaks, plot = FALSE)
  hist.data <- rbind(hist.data$mids, hist.data$counts)
  hist.data <- hist.data[ ,-c(1:2,(dim(hist.data)[2] - 2):dim(hist.data)[2])]
  result <- mean(hist.data[1, hist.data[2,] == min(hist.data[2,])])
  return(result)
}
.thresholdRanges <- function(x, percentage = 20){
  x <- x[!(x %in% NA)]
  result <- NULL
  result <- ((max(x) - min(x))/100) * percentage
  return(result)
}
.thresholdKmeans2 <- function(x, nClusters = 2, percentage = 30){
  x <- x[!(x %in% NA)]
  threshold <- kmeans(x = x, centers = nClusters)$centers
  if(dim(threshold)[1] == 2){
    threshold <- min(threshold) + (abs(diff(threshold)/100) * percentage)
    }
  return(threshold)
}
.thresholdDensityHist <- function(x, breaks = 100, strict = FALSE){
  x <- x[!(x %in% NA)]
  if(strict == TRUE){
    breaks <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
  }
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
.determineThresholds <- function(ch1 = NULL, ch2 = NULL, 
                                 algorithm = NULL, 
                                 breaks = 100, strict = TRUE,
                                 verbose = FALSE){
  if(is.null(ch1) == TRUE | is.null(ch2) == TRUE){ 
    stop("Data to determine the thresholds is not in correct format.\n")
  }
  if(is.null(algorithm) == TRUE){
    stop("algorithm has not been set to determine thresholds.\n")
  }
  ### REMOVE OUTLIERS
  if(verbose == TRUE){
    cat("Removing outliers from the data.\n")
  }
  amplitudes <- .removeOutliers(ch1, ch2, percentage = 1)
  ch1 <- amplitudes[,1]
  ch2 <- amplitudes[,2]
  ### RUN ALGORITHM
  if(tolower(algorithm) == "hist" | tolower(algorithm) == "histogram")
  {
    if(verbose == TRUE){
      cat("Setting threshold based on 'histogram'.\n")
    }
    result <- c(channel.1 = .thresholdHist(x = ch1, breaks = breaks, strict = strict), 
                channel.2 = .thresholdHist(x = ch2, breaks = breaks, strict = strict))
  }
  if(tolower(algorithm) == "ranges")
  {
    if(verbose == TRUE){
      cat("Setting threshold based on 'ranges'.\n")
    }
    result <- c(channel.1 = .thresholdRanges(x = ch1),
                channel.2 = .thresholdRanges(x = ch2))
  }
  if(tolower(algorithm) == "kmeans2")
  {
    if(verbose == TRUE){
      cat("Setting threshold based on 'kmeans'.\n")
    }
    result <- c(channel.1 = .thresholdKmeans2(x = ch1), 
                channel.2 = .thresholdKmeans2(x = ch2))
  }
  if(tolower(algorithm) == "densityhist")
  {
    if(verbose == TRUE){
      cat("Setting threshold based on 'densityhist'.\n")
    }
    result <- c(channel.1 = .thresholdDensityHist(x = ch1, breaks = breaks, strict = strict), 
                channel.2 = .thresholdDensityHist(x = ch2, breaks = breaks, strict = strict))
  }
 # result <- thresholdData(tData = tData, amplitude = result, type = 'threshold')
  return(result)
}
setThresholds <- function(data = NULL, algorithm = "densityhist", 
                          breaks = 100, strict = TRUE, 
                          type = "probe", verbose = TRUE){ 
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  }

  if(tolower(type) == "all"){
    channel.1 <- matrix(data = data@assayData$Ch1.Amplitude, ncol = 1)
    channel.2 <- matrix(data = data@assayData$Ch2.Amplitude, ncol = 1)
    result <- .determineThresholds(ch1 = channel.1, ch2 = channel.2, 
                         algorithm = algorithm, 
                         breaks = breaks, strict = strict,
                         verbose = verbose)
    data@phenoData$ch1['threshold',] <- result[1]
    data@phenoData$ch2['threshold',] <- result[2]
  } else if(tolower(type) == "probe"){
    probes <- unique(data@phenoData$sampleData['probe', ])
    for(i in 1:length(probes)){
      selection <- data@phenoData$sampleData['probe', ] %in% probes[i]
      channel.1 <- data@assayData$Ch1.Amplitude[,selection]
      channel.2 <- data@assayData$Ch2.Amplitude[,selection]
      channel.1 <- matrix(data = channel.1, ncol = 1)
      channel.2 <- matrix(data = channel.2, ncol = 1)
      result <- .determineThresholds(ch1 = channel.1, ch2 = channel.2, 
                                     algorithm = algorithm, 
                                     breaks = breaks, strict = strict,
                                     verbose = verbose)
      data@phenoData$ch1['threshold', selection] <- result[1]
      data@phenoData$ch2['threshold', selection] <- result[2]
    }
  }
  data <- .updateClusters(data) 
  return(data)
}
setThresholdsWell <- function(data = NULL, well = NULL, algorithm = "densityhist", 
                              breaks = 20, strict = TRUE, 
                              verbose = TRUE){ 
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  }
  if(is.null(well) != TRUE){
    if(class(well) == "numeric"){
      selection <- well
    } else if(class(well) == "character"){
      selection <- tolower(data@phenoData$sampleData['well',]) %in% tolower(well)
    }
  } else {
    stop("Well selection is not in the correct format.\n")
  }
    channel.1 <- data@assayData$Ch1.Amplitude[,selection]
    channel.2 <- data@assayData$Ch2.Amplitude[,selection]
    channel.1 <- matrix(data = channel.1, ncol = 1)
    channel.2 <- matrix(data = channel.2, ncol = 1)
    result <- .determineThresholds(ch1 = channel.1, ch2 = channel.2, 
                                     algorithm = algorithm, 
                                     breaks = breaks, strict = strict,
                                     verbose = verbose)
    data@phenoData$ch1['threshold', selection] <- result[1]
    data@phenoData$ch2['threshold', selection] <- result[2]
    data <- .updateClusters(data)
    return(data)
}
setThresholdsManual <- function(data = NULL, type = "probe",
                                verbose = TRUE){ 
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  }
  if(tolower(type) == "all"){
    channel.1 <- as.numeric(readline("Please enter numeric threshold value for channel 1: "))
    if(is.numeric(channel.1) != TRUE){
      stop("Threshold value has to be numeric. \n")
    }
    channel.2 <- as.numeric(readline("Please enter numeric threshold value for channel 2: "))
    if(is.numeric(channel.2) != TRUE){
      stop("Threshold value has to be numeric. \n")
    }

    data@phenoData$ch1['threshold',] <- channel.1
    data@phenoData$ch2['threshold',] <- channel.2
    
  } else if(tolower(type) == "probe"){
    probes <- unique(data@phenoData$sampleData['probe', ])
    for(i in 1:length(probes)){
      selection <- data@phenoData$sampleData['probe', ] %in% probes[i]
      message("Data for probe '", probes[i], "' can be entered.")
      cat("    type 'skip' if you do not wish to add a value for this probe.\n")
      
      channel.1 <- as.numeric(readline("Please enter numeric threshold value for channel 1: "))
      channel.2 <- as.numeric(readline("Please enter numeric threshold value for channel 2: "))
      if(is.numeric(channel.1) == TRUE & is.numeric(channel.2) == TRUE ){
        data@phenoData$ch1['threshold', selection] <- channel.1
        data@phenoData$ch2['threshold', selection] <- channel.2
        cat("   Threshold values have be updated. \n\n")
      }
    }
  }  
  data <- .updateClusters(data) 
  return(data)
}
