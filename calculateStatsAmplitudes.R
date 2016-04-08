calculateStatsAmplitudes <- function(x)
{
  col.names <- c("MeanAmplitudeofPositives","MeanAmplitudeofNegatives","MeanAmplitudeTotal")
  results <- matrix(0, nrow = 2, ncol = length(col.names),dimnames = list(NULL,col.names))
  results[1,colnames(results) == "MeanAmplitudeofPositives"] <- round(get.mean(x,cluster=c(2,3),1), digits = 2)
  results[2,colnames(results) == "MeanAmplitudeofPositives"] <- round(get.mean(x,cluster=c(4,3),2), digits = 2)
  results[1,colnames(results) == "MeanAmplitudeofNegatives"] <- round(get.mean(x,cluster=c(1,4),1), digits = 2)
  results[2,colnames(results) == "MeanAmplitudeofNegatives"] <- round(get.mean(x,cluster=c(1,2),2), digits = 2)
  results[1,colnames(results) == "MeanAmplitudeTotal"] <- round(get.mean(x,cluster=c(1,2,3,4),1), digits = 2)
  results[2,colnames(results) == "MeanAmplitudeTotal"] <- round(get.mean(x,cluster=c(1,2,3,4),2), digits = 2)
  return(results)
}