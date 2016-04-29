getThresholdsKmeans <- function(x, rm.outliers = TRUE)
{
  if(rm.outliers == TRUE)
  {
    x <- x[!x[,3] == 0,]
  }
  results <- kmeans(x[,1:2], 4)
  return(results)
}