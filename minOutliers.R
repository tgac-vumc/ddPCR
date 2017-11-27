minOutliers <- function(data = NULL, well = NULL, 
                        breaks = 250, strict = TRUE,
                        verbose = TRUE){
    if((class(data)[1] == "ddPCRdata") != TRUE){
      stop ("data structure is not in the correct format.\n\n")
    }
    if(verbose == TRUE){
      cat("Checking for outliers in each sample.\n")
    }
    if(is.null(well) == TRUE){
      wells <- 1:ncol(data@phenoData$sampleData)
    } else {
      wells <- well
    }
    
  for(z in 1:length(wells)){
    well <- wells[z]
    ### analyzing channel 1
    ch1 <- data@assayData$Ch1.Amplitude[ ,well]
    x <- ch1[!(ch1 %in% NA)]
    x <- x[x < data@phenoData$ch1['threshold', well]]
    
    if(strict == TRUE){
      breaks.ch1 <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
    } else {
      breaks.ch1 <- breaks
    }
    hist.data <- hist(x, breaks = breaks.ch1, plot = FALSE)
    x <- hist.data$counts
    result <- NULL
    for(i in 1:(length(x)-5)){
      if(x[i] > x[i+1] & x[i+1] <= x[i+2] &
         x[i+2] >= x[i+3] & x[i+3] >= x[i+4]) {
        result <- c(result, (i+1))
      }
    }
    
    
    if(is.null(result) != TRUE){
      result <- result[1]
      result.ch1 <- hist.data$mids[result]
      
      ### analyzing channel 2
      ch2 <- data@assayData$Ch2.Amplitude[ ,well]
      x <- ch2[!(ch2 %in% NA)]
      x <- x[x < data@phenoData$ch2['threshold', well]]
      
      if(strict == TRUE){
        breaks.ch2 <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
      } else {
        breaks.ch2 <- breaks
      }
      hist.data <- hist(x, breaks = breaks.ch2, plot = FALSE)
      x <- hist.data$counts
      result <- NULL
      for(i in 1:(length(x)-5)){
        if(x[i] > x[i+1] & x[i+1] <= x[i+2] &
           x[i+2] >= x[i+3] & x[i+3] >= x[i+4]) {
          result <- c(result, (i+1))
        }
      }
      
      if(is.null(result) != TRUE){
        result <- result[1]
        result.ch2 <- hist.data$mids[result]
        
        data@phenoData$ch1['minOutlier', well] <- result.ch1
        data@phenoData$ch2['minOutlier', well] <- result.ch2
      }
    }
  }
    data <- .updateClusters(data)
    return(data)
}
