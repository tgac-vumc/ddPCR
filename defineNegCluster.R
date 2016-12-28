getNegatives <- function(design){
  # get negatives from design file
}

defineNegCluster <- function(x, tData = NULL, design = NULL, well = NULL){

  # [ ] get negatives from design file
  if(is.null(design) == TRUE){
    getNegatives(design = design)
  }
  # [ ] get negatives from specific wells
  if(is.null(well) == TRUE){
    getPlateWells()
  }
  
  # [ ] combine negative data
  # [ ] remove outliers 5 %
  
  # [ ] calculate mean Sd for cluster 1
  
  return(tData)
}