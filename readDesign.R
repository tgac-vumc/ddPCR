readDesign <- function(path,pattern=NULL)
{
  design.file <- list.files(path=path,pattern=pattern)
  if(length(design.file) > 1){
    stop("multiple design files found. Unable to continue.\n"
    )
  }else{
    design <- read.table(file = file.path(path,design.file),header = TRUE,sep = "\t")
    return(design)
  }
}