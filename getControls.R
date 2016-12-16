getControls <- function(x, pos = c("positive", "pos control"), ntc = c("te buffer", "water"), neg = "")
{ # x = vector of sample names is input
  results <- rep("sample", length(x))
  x <- tolower(x)
  results[x %in% tolower(pos)] <- "pos"
  results[x %in% tolower(ntc)] <- "ntc"
  results[x %in% tolower(neg)] <- "neg"
  return(results)
}