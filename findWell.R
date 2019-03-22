.findWell <- function(files = NULL){
  if(is.null(files) == TRUE){
    stop("No files are given.\n")
  }
  result <- rep("NA",length(files))
  for(i in 1:length(files)){
    x <- strsplit(as.character(files[i]), split = "_")[[1]]
    result[i] <- x[x %in% .getPlateWells()]
  }
  return(result)
}