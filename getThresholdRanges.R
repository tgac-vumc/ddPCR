getThresholdRanges <- function(x)
{
  x <- as.numeric(x)
  result <- NULL
  result <- (max(x) - min(x)) / 2
  return(result)
}