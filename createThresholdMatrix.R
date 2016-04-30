createThresholdMatrix <- function()
{
  result <- matrix(data = NA, nrow = 7, ncol = 2, 
         dimnames = list(c("minOutlier","maxOutlier","threshold","thresholdMeanStDev","maxAmplitude","minRain","maxRain"), c("ch1","ch2")))
  return(result)
}
