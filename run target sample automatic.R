read.csv <- function(file){
  temp <- read.table(file = file,header = FALSE,sep = "\t",fill=TRUE)
  temp <- as.matrix(temp)
  nr.lines <- dim(temp)[1]
  nr.cols <- max(apply(X = temp,MARGIN = 1,FUN=countCharOccurrences))+1
  test.matrix <- matrix("NA",ncol=nr.cols,nrow=nr.lines)
  for (d in 1:(nr.lines)){
    x <-  unlist(strsplit(temp[d,1],split = ","))
    columns <- length(x)
    test.matrix[d,1:columns] <- x[1:columns]
  }
  return(test.matrix)
}

input.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\20151202 EGFR spike-in sm_2015-12-02-15-25"
error.file <- list.files(input.path, pattern = "Error.log",full.names = FALSE)
experiment.name <- gsub(pattern = "Error.log", replacement = "",x = error.file)
csv.file <- paste(experiment.name,".csv", sep="")
if(file.exists(file.path(input.path,csv.file)) == TRUE)
{
   experiment.data <-  read.csv(file = file.path(input.path,csv.file),header = TRUE, sep=",")
}
# - [ ]

targets <- unique(gsub(pattern = " wt", replacement = "", x=experiment.data$Target))
wells <- unique(experiment.data$Well)

file.exists(
file.path(input.path,paste(experiment.name,"_",wells,"_Amplitude.csv",sep=""))
)

if(toupper(target) == "L858R")
{
  control.sample <- "H1975"
}
if(toupper(target) == "T790M")
{
  control.sample <- "H1975"
}

get.controls <- function(x, pos=c("positive","pos control"), ntc=c("te buffer","water"), neg="")
{
  results <- rep("sample",length(x))
  x <- tolower(x)
    results[x %in% tolower(pos)] <- "pos"
    results[x %in% tolower(ntc)] <- "ntc"
    results[x %in% tolower(neg)] <- "neg"
    return(results)
}

# loop for targets
# if all files exist
# sub text for labels  = Target 1 vs target 2 = ie L858R & L858R wt

get.experiments <- function(path)
{
  file <- list.files(path, pattern = "Error.log",full.names = TRUE)
  file <- gsub(pattern = "Error.log", replacement = "",x = file)
  experiment.name <- basename(file)
  file <- paste(file,".csv", sep="")
  if(file.exists(file) == TRUE)
  {
    data <-  read.table(file = file,header = TRUE, sep=",",row.names = NULL)
  } else { break 
    }
  targets <- unique(gsub(pattern = " wt", replacement = "", x=data$Target))
  grep(data$Target,pattern = targets[i])
  
  result <- list()
  for(i in 1:length(targets))
  {
      result[[i]] <- data[grep(data$Target,pattern = targets[i]),]
  }
  names(result) <- as.character(targets)
  return(result)
}

