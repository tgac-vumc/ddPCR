dropletCountText <- function(x)
{
  results <- NULL
  results <- list(clusters = c(cluster.1 = sum(x == 1), cluster.2 = sum(x == 2), cluster.3 = sum(x == 3), cluster.4 = sum(x == 4)))
  results <- paste("Ch1-Ch2-:", results$clusters[1],
                   "   Ch1+Ch2-:", results$clusters[2],
                   "   Ch1+Ch2+:", results$clusters[3],
                   "   Ch1-Ch2+:", results$clusters[4], sep="")
  return(results)
}