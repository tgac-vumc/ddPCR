createThresholdMatrix <- function(){
  result <- matrix(data = NA, nrow = 7, ncol = 2, 
         dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minAmplitude","minRain","maxRain"), c("ch1","ch2")))
  return(result)
}
addThresholdData <- function(x = NULL, add, type = ''){
  
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