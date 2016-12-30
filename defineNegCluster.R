removeOutliers <- function(x, cutoff = 5){
  percentage <- round((nrow(x)/100) * (cutoff/2))
  selection <- c(1:percentage, (nrow(x)-percentage):nrow(x))
  selection <- !(1:nrow(x)  %in% selection)
  x <- x[order(x[,1]),]
  x <- x[selection,]
  selection <- c(1:percentage, (nrow(x)-percentage):nrow(x))
  selection <- !(1:nrow(x)  %in% selection)
  x <- x[order(x[,2]),]
  x <- x[selection,]
  return(x)
}
getNegatives <- function(design = NULL, verbose = TRUE){
  if(is.null(design) == TRUE){
      stop("No design file is given.\n")
  } else {
    negatives <- design$Type == "neg"
    if(TRUE %in% negatives == TRUE){
      design <- design[negatives,]
    } else {
      stop("No negative control was found.\n")
    }
  } 
  return(design)
}
defineNegCluster <- function(x, tData = NULL, design = NULL, well = NULL){

  # [ ] get negatives from design file
  if(is.null(design) == FALSE){
    getNegatives(design = design)
  }
  # [ ] get negatives from specific wells
  if(is.null(well) == FALSE){
    getPlateWells()
  }
  
  # [ ] combine negative data
  # [ ] remove outliers 5 %
  
  # [ ] calculate mean Sd for cluster 1
  
  return(tData)
}

# [ ] function: highest density
# [ ] function: get negatives from design
# [ ] create pipeline for neg analysis
# [ ] 
# [ ] 