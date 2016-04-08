defineClusters <- function(x, thresholds)
{
  if (length(thresholds) != 2) {
    stop("thresholds must have a length of 2.\n")
  }
  results <- rep(NA,dim(x)[1])
  results[x[,1] < thresholds[1] & x[,2] < thresholds[2]] <- 1 # ch1-ch2- : cluster 1
  results[x[,1] > thresholds[1] & x[,2] < thresholds[2]] <- 2 # ch1+ch2- : cluster 2
  results[x[,1] > thresholds[1] & x[,2] > thresholds[2]] <- 3 # ch1+ch2+ : cluster 3
  results[x[,1] < thresholds[1] & x[,2] > thresholds[2]] <- 4 # ch1-ch2+ : cluster 4
  x[,3] <- results
  return(x)
}