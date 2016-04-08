getTargets <- function(path)
{ # will use .csv file of the experiment located in path
  file <- list.files(path, pattern = "Error.log",full.names = TRUE)
  file <- gsub(pattern = "Error.log", replacement = "",x = file)
  experiment.name <- basename(file)
  file <- paste(file,".csv", sep="")
  if(file.exists(file) == TRUE)
  {
    data <-  read.csv(file = file,header = TRUE,sep = ",",check.names = FALSE, row.names = NULL)
    if(colnames(data)[1] == "row.names")
    {
      colnames(data) <- colnames(data)[2:dim(data)[2]]
      data <- data[,1:(dim(data)[2]-1)]
    }
  } else { break 
  }
  targets <- unique(mgsub(pattern = c("_WT","_wt", " wt"," WT"), replacement = c("","","",""), x=data$Target))
  result <- list()
  for(i in 1:length(targets))
  {
    result[[i]] <- data[grep(data$Target,pattern = targets[i]),]
  }
  targets <- gsub(pattern = " ", replacement = "_", x=targets)
  names(result) <- targets
  return(result)
}