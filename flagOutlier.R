flagOutlier <- function(data = NULL, well = NULL){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  } 
  if(is.null(well) == TRUE){
    cat("No well location was given to plot.\n")
  } 
  if(class(well) == "integer"){
    well <- as.numeric(well)
  }
  if (class(well) == "character"){
    sample <- match(tolower(well), table = tolower(data@phenoData$sampleData['well',]))
    if(is.na(sample) == TRUE){
      cat("No correct sample well has been given.\n")
    }
  } else if(class(well) == "numeric"){
    if((well > ncol(data@phenoData$sampleData)) == TRUE){
      cat("No correct sample well has been given.\n")
    }else {
      sample <- well
    }
  } 
  data@phenoData$sampleData['probe', well] <- "outlier"
  return(data)
}
