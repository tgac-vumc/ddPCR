convertCopiesToNanogram <- function(x)
{
  x <- as.numeric(x) / 15.152
  return(x)
}