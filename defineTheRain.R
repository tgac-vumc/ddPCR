defineTheRain <- function(x, stdev = c(3,3), tData = NULL)
{
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