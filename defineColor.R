defineColor <- function(x, density = NULL)
{
  ddpcr.colors <- paste(c("#000000", "#FF6600", "#00CC00", "#0033FF", "#FF0000" , "#FF00CC"), sep = "")
  x <- as.character(x)
  x <- mgsub(pattern = c("1","3","4","2", "5", "0"), replacement = ddpcr.colors, x = x)
  if(class(density) != "NULL")
  {
    x <- paste(x, as.character(density), sep = "") 
  }
  
  return(x)
}