minOutliers <- function(data = NULL, well = NULL, breaks = 250, strict = FALSE){
    if((class(data)[1] == "ddPCRdata") != TRUE){
      stop ("data structure is not in the correct format.\n\n")
    }
    if(is.null(well) == TRUE){
      wells <- ncol(data@phenoData$sampleData)
    } else {
      wells <- well
    }
    
  for (i in length(wells)){
    well <- wells[i]
    ch1 <- data@assayData$Ch1.Amplitude[ ,well]
    x <- ch1[!(ch1 %in% NA)]
    if(strict == TRUE){
      breaks <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
    }
    hist.data <- hist(x, breaks = breaks, plot = FALSE)
    x <- hist.data$counts
    result <- NULL
    for(i in length(x):6){
      if(x[i-1] < x[i-2] & x[i-2] < x[i-3] &
         x[i-3] >= x[i-4] & x[i-4] >= x[i-5]) {
        result <- c(result, (i-3))
      }
    }
    result <- result[length(result)]
    result.ch1 <- hist.data$mids[result]
    
    ch2 <- data@assayData$Ch2.Amplitude[ ,well]
    x <- ch2[!(ch2 %in% NA)]
    if(strict == TRUE){
      breaks <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
    }
    hist.data <- hist(x, breaks = breaks, plot = FALSE)
    x <- hist.data$counts
    result <- NULL
    for(i in length(x):6){
      if(x[i-1] < x[i-2] & x[i-2] < x[i-3] &
         x[i-3] >= x[i-4] & x[i-4] >= x[i-5] ) {
        result <- c(result, (i-3))
      }
    }
    result <- result[length(result)]
    result.ch2 <- hist.data$mids[result]
    
    data@phenoData$ch1['minOutlier', well] <- result.ch1
    data@phenoData$ch2['minOutlier', well] <- result.ch2
  }
    data <- .updateClusters(data)
    return(data)
  }