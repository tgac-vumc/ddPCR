calculateStatsDroplets <- function(x)
{ # input = amplitide data with defined clusters
  col.names <- c("Positives","Negatives","Ch1-Ch2-","Ch1+Ch2-","Ch1+Ch2+","Ch1-Ch2+","AcceptedDroplets")
  results <- matrix(0, nrow=2,ncol=length(col.names),dimnames = list(NULL,col.names))
  results[1,colnames(results) == "Positives"] <- droplet.count(x, c(2,3)) #Positives
  results[2,colnames(results) == "Positives"] <- droplet.count(x, c(3,4)) #Positives
  results[1,colnames(results) == "Negatives"] <- droplet.count(x, c(1,4))# Negatives
  results[2,colnames(results) == "Negatives"] <- droplet.count(x, c(1,2)) # Negatives
  results[1:2,colnames(results) == "Ch1-Ch2-"] <- droplet.count(x, 1)
  results[1:2,colnames(results) == "Ch1+Ch2-"] <- droplet.count(x, 2)
  results[1:2,colnames(results) == "Ch1+Ch2+"] <- droplet.count(x, 3)
  results[1:2,colnames(results) == "Ch1-Ch2+"] <- droplet.count(x, 4)
  results[1:2,colnames(results) == "AcceptedDroplets"] <- droplet.count(x, c(1,2,3,4)) #AcceptedDroplets
  return(results)
}