combineSamples <- function(path,files)
{
  combined.data <- NULL
  for(i in 1:length(files))
  { 
    sample.data <- read.table(file=file.path(path,files[i]),header = TRUE,sep = ",")
    combined.data <- rbind(combined.data,sample.data )
  }
  return(combined.data)
}