thresholdData <- function(tData = NULL, amplitude = NULL, type = NULL) # tData = thresholdData
{
  if(class(type) != "NULL")
  {
    thresholds <- c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minRain","maxRain")
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