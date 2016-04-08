getThresholdsKmeans <- function(x)
{
  results <- kmeans(x[,1:2], 4)
  return(results)
}