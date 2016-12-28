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
