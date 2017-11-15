getMaxAmplitude <- function(x, tData = NULL)
{
  result <- c(round(max(x[,1]) + 100), round(max(x[,2]) + 100))
  result <- thresholdData(tData = tData, amplitude = result, type = 'maxAmplitude')
  return(result)
}
getMinAmplitude <- function(x, tData = NULL)
{
  result <- c(round(min(x[,1]) + 100), round(min(x[,2]) + 100))
  result <- thresholdData(tData = tData, amplitude = result, type = 'minAmplitude')
  return(result)
}