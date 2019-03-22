changeColors <- function(data, colors = NULL){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n\n")
  } 
  if(is.null(colors) == TRUE){
    stop ("Colors are not in the correct format.\n\n")
  } else if(length(colors) < 6 | length(colors) > 6){
    stop ("Colors must be a length of 6.\n\n")
  }
  if(sum(substr(colors, 1, 1) %in% "#") < 6){
    stop ("Colors must be in hex format.\n\n")
  }
  
  data@experimentData$Cluster.Color <- colors
  data <- .updateColors(data)  
  return(data)
}