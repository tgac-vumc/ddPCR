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
.thresholdHist <- function(x, breaks = 20, strict = FALSE){
  x <- x[!(x %in% NA)]
  cutoff <- min(x) + (abs(max(x) - min(x)) / 2)
  if(strict == TRUE){
    breaks <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
  }
  hist.data <- hist(x, breaks = breaks, plot = FALSE)
  hist.data <- rbind(mids = hist.data$mids, counts = hist.data$counts)
  hist.data <- hist.data[ ,hist.data[1,] < cutoff]
  
  hist.data <- hist.data[,(grep(pattern = max(hist.data[2,]), x = hist.data[2,])):dim(hist.data)[2]]
  result <- mean(hist.data[1, hist.data[2,] == min(hist.data[2,])])
  return(result)
}
.thresholdRanges <- function(x, percentage = 33){
  x <- x[!(x %in% NA)]
  result <- NULL
  result <- ((max(x) - min(x))/100) * percentage
  result <- result + min(x)
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
.thresholdDensityHist <- function(x, breaks = 100, strict = FALSE, verbose = TRUE){
  x <- x[!(x %in% NA)]
  # remove high positive droplets for analysis
  cutoff <- min(x) + (abs(max(x) - min(x)) / 2)
  if(strict == TRUE){
    breaks <- .makeBreaks(min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE), breaks = breaks)
  }
  hist.data <- hist(x, breaks = breaks, plot = FALSE)
  x <- hist.data$counts
  x <- c(x, x[length(x)], x[length(x)])
  result <- NULL
  for(i in 1:(length(x)-3)){
    if(x[i+1] > x[i+2] & x[i+2] > x[i+3] &
       x[i+3] <= x[i+4]
    ) {
      result <- c(result, i+3)
    }
  }
  result <- hist.data$breaks[result]
  result <- result[result < cutoff]
  if(length(result) > 1){
    result <- max(result)
  }
  return(result)
}
.thresholdmeanSd <- function(x, stdev = 3, breaks = 15){
  x <- x[!(x %in% NA)]
  threshold <- .thresholdHist(x = x, breaks = breaks)
  x <- x[x < threshold]
  refined.threshold <- mean(x) + (stdev * sd(x))
  if(refined.threshold > threshold){
    refined.threshold <- threshold
  } 
  return(refined.threshold)
}
.thresholdMeanDensityMirror <- function(x, breaks = 100, strict = FALSE){
  x <- x[!(x %in% NA)]
  if(strict == TRUE){
    breaks <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
  }
  hist.data <- hist(x, breaks = breaks, plot = FALSE)
  highest.desity <- hist.data$mids[hist.data$counts == max(hist.data$counts)]
  result <- diff(x = c(min(x), highest.desity)) + highest.desity
  
  return(result)
}
.determineThresholds <- function(ch1 = NULL, ch2 = NULL, 
                                 algorithm = NULL, 
                                 breaks = 100, strict = TRUE,
                                 rm.percentage = NULL,
                                 verbose = FALSE){
  if(is.null(ch1) == TRUE | is.null(ch2) == TRUE){ 
    stop("Data to determine the thresholds is not in correct format.\n")
  }
  if(is.null(algorithm) == TRUE){
    stop("algorithm has not been set to determine thresholds.\n")
  }
  
  ### removing percentage
  if(is.numeric(rm.percentage) == TRUE){
    ch1 <- .removePercentage(ch1, percentage = rm.percentage)
    ch2 <- .removePercentage(ch2, percentage = rm.percentage)
  }
  
  ### RUN ALGORITHM
  if(tolower(algorithm) == "hist" | tolower(algorithm) == "histogram")
  {
    result <- c(channel.1 = .thresholdHist(x = ch1, breaks = breaks, strict = strict), 
                channel.2 = .thresholdHist(x = ch2, breaks = breaks, strict = strict))
  }
  if(tolower(algorithm) == "ranges")
  {
    result <- c(channel.1 = .thresholdRanges(x = ch1),
                channel.2 = .thresholdRanges(x = ch2))
  }
  if(tolower(algorithm) == "kmeans2")
  {
    result <- c(channel.1 = .thresholdKmeans2(x = ch1), 
                channel.2 = .thresholdKmeans2(x = ch2))
  }
  if(tolower(algorithm) == "densityhist")
  {
    result <- c(channel.1 = .thresholdDensityHist(x = ch1, breaks = breaks, strict = strict), 
                channel.2 = .thresholdDensityHist(x = ch2, breaks = breaks, strict = strict))
  }
  if(tolower(algorithm) == "meansd")
  {
    result <- c(channel.1 = .thresholdmeanSd(x = ch1), 
                channel.2 = .thresholdmeanSd(x = ch2))
  }
  if(tolower(algorithm) == "mirror")
  {
    result <- c(channel.1 = .thresholdMeanDensityMirror(x = ch1), 
                channel.2 = .thresholdMeanDensityMirror(x = ch2))
  }
  return(result)
}
setThresholds <- function(data = NULL, algorithm = "densityhist", 
                          breaks = 20, strict = TRUE, 
                          type = "probe",
                          rm.percentage = NULL,
                          verbose = TRUE){ 
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  }
  
  
  ### analyse all samples together
  if(tolower(type) == "all"){
    ### get data from matrix
    channel.1 <- matrix(data = data@assayData$Ch1.Amplitude, ncol = 1)
    channel.2 <- matrix(data = data@assayData$Ch2.Amplitude, ncol = 1)
    
    result <- .determineThresholds(ch1 = channel.1, ch2 = channel.2, 
                         algorithm = algorithm, 
                         breaks = breaks, strict = strict,
                         rm.percentage = rm.percentage,
                         verbose = verbose)
    
    data@phenoData$ch1['threshold',] <- result[1]
    data@phenoData$ch2['threshold',] <- result[2]
  }
  
  ### analyse samples per probe
  if(tolower(type) == "probe"){
    
    if(verbose == TRUE){
      cat("Setting threshold based on '", algorithm, "' for each probeset.\n", sep="")
    }
    ### find all the unique probes
    probes <- unique(data@phenoData$sampleData['probe', ])
    for(i in 1:length(probes)){
      selection <- data@phenoData$sampleData['probe', ] %in% probes[i]
      
      channel.1 <- data@assayData$Ch1.Amplitude[,selection]
      channel.2 <- data@assayData$Ch2.Amplitude[,selection]
      channel.1 <- matrix(data = channel.1, ncol = 1)
      channel.2 <- matrix(data = channel.2, ncol = 1)

      ### run analysis 
      result <- .determineThresholds(ch1 = channel.1, ch2 = channel.2, 
                                     algorithm = algorithm, 
                                     breaks = breaks, strict = strict,
                                     rm.percentage = rm.percentage,
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
                              rm.percentage = NULL,
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
                                     rm.percentage = rm.percentage,
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
