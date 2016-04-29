poissonCorrection <- function(posCount, count, iterations=10)
{
  if(posCount == 0)
  {
    return(0)
  } else 
  {
    additional.droplets <- rep(NA, iterations)
    for(i in 1:iterations)
    {
      result <- rep(NA, count)
      droplets <- 1:count
      result[1:posCount] <- sample(x = droplets, size = posCount, replace = TRUE)
      counter <- posCount + 1
      for(z in (posCount + 1):count)
      {
        if((length(unique(result))-1) >= posCount)
        {
          break
        }
        result[z] <- sample(x = droplets,size = 1, replace = TRUE)
        z <- z+1
      }
      additional.droplets[i] <- z - posCount
    }
    additional.droplets <- mean(c((posCount + min(additional.droplets)), posCount + max(additional.droplets)))
    return(additional.droplets)
  }
}