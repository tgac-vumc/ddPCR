defineColor <- function(x, density = NULL)
{
  result <- rep(NA, length(x))
  result[grep(pattern = 0, x = x)] <- "#FF0000"
  result[grep(pattern = 1, x = x)] <- "#000000"
  result[grep(pattern = 2, x = x)] <- "#0033FF" 
  result[grep(pattern = 3, x = x)] <- "#FF6600"
  result[grep(pattern = 4, x = x)] <- "#00CC00"
  result[grep(pattern = 5, x = x)] <- "#FF00CC"
  if(class(density) != "NULL")
  {
    result <- paste(result, as.character(density), sep = "") 
  }
  return(result)
}