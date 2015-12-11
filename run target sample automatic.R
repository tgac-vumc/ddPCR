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


if(toupper(target) == "L858R")
{
  control.sample <- "H1975"
}
if(toupper(target) == "T790M")
{
  control.sample <- "H1975"
}



# loop for targets
# if all files exist
# sub text for labels  = Target 1 vs target 2 = ie L858R & L858R wt
# - [ ] statistics: multiple small functions

get.statistics <- function(x,sample=NULL,input.file=NULL,target=NULL,breakpoints=NULL, interations=100)
{

  results <- matrix(NA, nrow = 2, ncol = length(col.names),dimnames = list(NULL,col.names))
  if(class(input.file) != "NULL"){ results[1:2,colnames(results) == "Well"] <- get.well(input.file)}
  if(class(sample) != "NULL"){results[1:2,colnames(results) == "Sample"] <- sample} 
  results[1:2,colnames(results) == "TargetType"] <- c("Ch1","Ch2")
  if(class(target) != "NULL"){ results[1:2,colnames(results) == "Target"] <- sample}
  

  
  if(length(breakpoints) == 2)
  {
    results[1:2,colnames(results) == "Threshold"] <- breakpoints
  }

  results[1:2,colnames(results) == "ConcentrationPer1ul"] <- as.numeric(results[1:2,colnames(results) == "CopiesPer1ul"])/15.152
  


  
get.statistics.quality <- function(x,y)
{
  col.names <- c("Status")
  results <- matrix(NA, nrow = 2, ncol = length(col.names),dimnames = list(NULL,col.names))
  results[1:2,colnames(results) == "Status"] <- "OK"
  
  channel.1 <- round(calc.copies(posCount = as.numeric(results[1,colnames(results) == "Positives"]),count = as.numeric(results[1,colnames(results) == "AcceptedDroplets"])), digits = 1)
  channel.2 <- round(calc.copies(posCount = as.numeric(results[2,colnames(results) == "Positives"]),count = as.numeric(results[2,colnames(results) == "AcceptedDroplets"])), digits = 1)
  
  
  
  if(as.numeric(results[1:2,colnames(results) == "AcceptedDroplets"]) < 10000){results[1:2,colnames(results) == "Status"] <- "CHECK"}
  if(as.numeric(results[1,colnames(results) == "CopiesPer1ul"]) == 0){results[1:2,colnames(results) == "Status"] <- "CHECK"}
  if(as.numeric(results[2,colnames(results) == "CopiesPer1ul"]) == 0){results[1:2,colnames(results) == "Status"] <- "CHECK"}
  return(results)
}

}

