.minOutliers <- function(data = NULL, breaks = 300, strict = TRUE)
{
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n\n")
  } 
  # function needs to be converted for S4 structure.
  
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

