defineColor <- function(x,density=NULL)
{
  ddpcr.colors <- paste(c("#000000","#FF6600","#00CC00","#0033FF"), sep="")
  x <- as.character(x)
  x <- mgsub(pattern = c("1","3","4","2"),replacement = ddpcr.colors,x=x)
  if(class(density) != "NULL")
  {
    x <- paste(x, as.character(density),sep="") 
  }
  
  return(x)
}