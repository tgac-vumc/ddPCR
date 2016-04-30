getMaxAmplitude <- function(x)
{
  result <- c(round(max(x[,1]) + 100), round(max(x[,2]) + 100))
  result <- matrix(data = result, nrow = 1, ncol = 2, dimnames = list(c("maxAmplitude"), c("ch1","ch2")))
  return(result)
}