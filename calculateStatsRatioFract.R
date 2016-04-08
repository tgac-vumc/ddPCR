calculateStatsRatioFract <- function(x)
{ # input = output get.statistics.copies()
  if ("CopiesPer1ul" %in% colnames(x) == TRUE)
  {
    col.names <- c("Ratio","FractionalAbundance")
    results <- matrix(0, nrow = 2, ncol = length(col.names),dimnames = list(NULL,col.names))
    results[1:2,colnames(results) == "Ratio"] <- x[1,colnames(x) == "CopiesPer1ul"] / x[2,colnames(x) == "CopiesPer1ul"]
    results[2,colnames(results) == "Ratio"] <- x[2,colnames(x) == "CopiesPer1ul"] / x[1,colnames(x) == "CopiesPer1ul"]
    results[1,colnames(results) == "FractionalAbundance"] <- x[1,colnames(x) == "CopiesPer1ul"] / (x[1,colnames(x) == "CopiesPer1ul"]+x[2,colnames(x) == "CopiesPer1ul"])*100
    results[2,colnames(results) == "FractionalAbundance"] <- x[2,colnames(x) == "CopiesPer1ul"] / (x[1,colnames(x) == "CopiesPer1ul"]+x[2,colnames(x) == "CopiesPer1ul"])*100
    return(results)
  }
}