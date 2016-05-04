getMaxAmplitude <- function(x, tData = NULL)
{
  result <- c(round(max(x[,1]) + 100), round(max(x[,2]) + 100))
  result <- thresholdData(tData = tData, amplitude = result, type = 'maxAmplitude')
  return(result)
}