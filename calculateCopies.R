calculateCopies <-  function(posCount, count, vDroplet=0.91, volume=1)
{ # concentration in copies / user defined volume
  negCount <- count-posCount
  result <- ((-log(negCount/count)/vDroplet)) * 1000 * volume
  return(result)
}