getMaxAmplitude <- function(x)
{
  results <- c(Ch1.max = round(max(x[,1]) + 100), Ch2.max = round(max(x[,2]) + 100))

  return(results)
}