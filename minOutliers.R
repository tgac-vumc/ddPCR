minOutliers <- function(data = NULL, breaks = 250, strict = TRUE){
  answer <- readline(prompt <- "function is not working correctly. Type 'YES' to continue : ")
  if(tolower(answer) == "y" | tolower(answer) == "yes" ){
    cat("Ignoring warning..\n")
  } else (stop("stopped manually.\n"))

  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n\n")
  }
  for (i in 1:ncol(data@phenoData$sampleData)){
    x <- data@assayData$Ch1.Amplitude[ ,i]
    x <- x[!(x %in% NA)]
    if(strict == TRUE){
      breaks.ch1 <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
    } else { 
      breaks.ch1 <- breaks
      }
    hist.data <- hist(x, breaks=breaks.ch1, plot=FALSE)
    firstBigClusterCh1 <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
    if(TRUE %in% (hist.data$counts[1:firstBigClusterCh1] %in% 0 == TRUE))
      {
      closestZero <- max((1:firstBigClusterCh1)[hist.data$counts[1:firstBigClusterCh1] == 0])
      minOutlierCh1 <- hist.data$mids[closestZero]
    } else  { minOutlierCh1 <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) 
              if(minOutlierCh1 < 0){ minOutlierCh1 <- 0 }
            }
    
    x <- data@assayData$Ch2.Amplitude[ ,i]
    x <- x[!(x %in% NA)]
    if(strict == TRUE){
      breaks.ch2 <- .makeBreaks(min = min(x), max = max(x), breaks = breaks)
    } else { 
      breaks.ch2 <- breaks
    }
    hist.data <- hist(x, breaks=breaks.ch2, plot=FALSE)
    firstBigClusterCh2 <- grep(pattern = max(hist.data$counts), x = hist.data$counts)
    if(TRUE %in% (hist.data$counts[1:firstBigClusterCh2] %in% 0 == TRUE))
    {
      closestZero <- max((1:firstBigClusterCh2)[hist.data$counts[1:firstBigClusterCh2] == 0])
      minOutlierCh2 <- hist.data$mids[closestZero]
    } else { minOutlierCh2 <- hist.data$breaks[1] - diff(hist.data$breaks[1:2]) 
            if(minOutlierCh2 < 0){ minOutlierCh2 <- 0 }
            }
    
    data@phenoData$ch1['minOutlier', i] <- minOutlierCh1
    data@phenoData$ch2['minOutlier', i] <- minOutlierCh2
  }
  data <- .updateClusters(data)
  return(data)
}
