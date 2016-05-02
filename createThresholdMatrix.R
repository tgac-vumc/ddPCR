createThresholdMatrix <- function()
{
  result <- matrix(data = NA, nrow = 7, ncol = 2, 
         dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minRain","maxRain"), c("ch1","ch2")))
  return(result)
}
addThresholdData <- function(x = NULL, add, type = '')
{
  
  if (class(result != "matrix"))
  {
    result <- matrix(data = NA, nrow = 7, ncol = 2, 
                     dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minRain","maxRain"), c("ch1","ch2")))
    
  }
  dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minRain","maxRain"), c("ch1","ch2")))

}
# [ ] - create new function that will create matrix or add settings to it.
# some fixed settings
# able to add own settings as well
# 