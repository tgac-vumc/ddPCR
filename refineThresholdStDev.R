refineThresholdStDev <- function(x, stdev=3, thresholds)
{
  cluster1 <- colSums(calculaeMeanSdCluster(x, cluster = 1, stdev = stdev))
  cluster2 <- colSums(calculateMeanSdCluster(x, cluster = 2, stdev = stdev))
  cluster4 <- colSums(calculateMeanSdCluster(x, cluster = 4, stdev = stdev))
  refined.channel2 <- max(c(cluster1[2], cluster2[2]), na.rm = TRUE)
  if(refined.channel2 < thresholds[2])
  {
    thresholds[2] <- mean(c(refined.channel2, thresholds[2]))
  } else {cat("threshold for channel 2 could not be defined further.\n")}
  refined.channel1 <- max(c(cluster1[1], cluster4[1]), na.rm = TRUE)
  if(refined.channel1 < thresholds[1])
  {
    thresholds[1] <- mean(c(refined.channel1, thresholds[1]))
  } else {cat("threshold for channel 1 could not be defined further.\n")}
  return(thresholds)
}