.updateClusters <- function(data = NULL){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n\n")
  } 
  
  for(i in 1:ncol(data@phenoData$sampleData)){
    if(is.na(data@phenoData$ch1['threshold', i]) != TRUE){
      threshold <- c(channel.1 = data@phenoData$ch1['threshold', i],
                     channel.2 = data@phenoData$ch2['threshold', i])

      x <- cbind(channel.1 = data@assayData$Ch1.Amplitude[ ,i], 
                 channel.2 = data@assayData$Ch2.Amplitude[ ,i])
      
      results <- data@assayData$Cluster[, i]
      results[x[,1] < threshold[1] & x[,2] < threshold[2]] <- 1 # ch1-ch2- : cluster 1
      results[x[,1] > threshold[1] & x[,2] < threshold[2]] <- 2 # ch1+ch2- : cluster 2
      results[x[,1] > threshold[1] & x[,2] > threshold[2]] <- 3 # ch1+ch2+ : cluster 3
      results[x[,1] < threshold[1] & x[,2] > threshold[2]] <- 4 # ch1-ch2+ : cluster 4
      
      data@assayData$Cluster[, i] <- results
    }
    if(is.na(data@phenoData$ch1['minOutlier', i]) != TRUE){
      
      minOutlier <- c(channel.1 = data@phenoData$ch1['minOutlier', i],
                     channel.2 = data@phenoData$ch2['minOutlier', i])
      threshold <- c(channel.1 = data@phenoData$ch1['threshold', i],
                     channel.2 = data@phenoData$ch2['threshold', i])
      
      x <- cbind(channel.1 = data@assayData$Ch1.Amplitude[ ,i], 
                 channel.2 = data@assayData$Ch2.Amplitude[ ,i])
      
      results <- data@assayData$Cluster[, i]
      results[x[,1] < minOutlier[1] & x[,2] < threshold[2]] <- 0 
      results[x[,2] < minOutlier[2] & x[,1] < threshold[1]] <- 0 
      
      data@assayData$Cluster[, i] <- results
    }
    if(is.na(data@phenoData$ch1['maxOutlier', i]) != TRUE){
      
      maxOutlier <- c(channel.1 = data@phenoData$ch1['maxOutlier', i],
                      channel.2 = data@phenoData$ch2['maxOutlier', i])
      
      x <- cbind(channel.1 = data@assayData$Ch1.Amplitude[ ,i], 
                 channel.2 = data@assayData$Ch2.Amplitude[ ,i])
      
      results <- data@assayData$Cluster[, i]
      results[x[,1] > maxOutlier[1]] <- 0
      results[x[,2] > maxOutlier[2]] <- 0 
      
      data@assayData$Cluster[, i] <- results
    }
  }
  data <- .updateColors(data)
  return(data)
}
