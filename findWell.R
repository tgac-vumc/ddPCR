findWell <- function(input.file)
{
  x <- strsplit(as.character(input.file), split = "_")[[1]]
  x <- x[x %in% getPlateWells()]
  return(x)
}