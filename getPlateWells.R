.getPlateWells <- function(prefix=NULL, suffix=NULL){
  result <- rep(NA, 96)
  
  numbers <- c("01","02","03","04","05","06","07","08","09","10","11","12")
  start <- 1
  for(i in 1:8)
  {
    result[start:(start + 11)] <- paste(LETTERS[i], numbers, sep="")
    start <- start + 12
  }
  if(class(prefix) != "NULL")
  {
    result <- paste(prefix, result, sep="")
  }
  if(class(suffix) != "NULL")
  {
    result <- paste(result, suffix, sep="")
  }
  return(result)
}