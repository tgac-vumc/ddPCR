.updateColors <- function(data, density = NULL){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n\n")
  } 
  
  Cluster <- data@assayData$Cluster
  Color <- data@assayData$Color
  
  for(i in 1:length(data@experimentData$Cluster.Number)){
    selection <- Cluster %in% data@experimentData$Cluster.Number[i]
    Color[selection] <- data@experimentData$Cluster.Color[i]
  }
  Color <- apply(X = Color, FUN = paste0, MARGIN = 2, ... = density)
  Color[Color %in% paste0(NA, density)] <- NA
  data@assayData$Color <- Color
  
  return(data)
}