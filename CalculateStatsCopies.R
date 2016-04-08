CalculateStatsCopies <- function(x, interations = 100)
{ # input = amplitide data with defined clusters
  col.names <- c("CopiesPer1ul","CopiesPer20ulWell","PoissonCopiesPer1ul")
  results <- matrix(0, nrow = 2, ncol = length(col.names),dimnames = list(NULL,col.names))
  channel.1 <- droplet.count(x, c(2,3)) #Positives channel 1
  channel.2 <- droplet.count(x, c(3,4)) #Positives channel 1
  total <- droplet.count(x, c(1,2,3,4)) #AcceptedDroplets
  results[1,colnames(results) == "CopiesPer1ul"] <- round(calc.copies(posCount = channel.1, count = total), digits = 1)
  results[2,colnames(results) == "CopiesPer1ul"] <- round(calc.copies(posCount = channel.2, count = total), digits = 1)
  results[1,colnames(results) == "PoissonCopiesPer1ul"] <- round(calc.copies(posCount = poisson.correction(count = total, posCount = channel.1, iterations = interations), count = total), digits = 1)
  results[2,colnames(results) == "PoissonCopiesPer1ul"] <- round(calc.copies(posCount = poisson.correction(count = total, posCount = channel.2, iterations = interations), count = total), digits = 1)
  results[1:2,colnames(results) == "CopiesPer20ulWell"] <- results[1:2,colnames(results) == "CopiesPer1ul"]*20
  return(results)
}