samples <- function(data){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  } 
  print(data@phenoData$sampleData)
}