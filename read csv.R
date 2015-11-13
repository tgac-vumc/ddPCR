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