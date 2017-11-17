checkPackages <- function(x){
  for (i in 1:length(x)){
    package <- as.character(x[i])
    if( package %in% installed.packages()[,1] == FALSE ){
      install.packages(package, lib = "N:/Documenten/R/win-library/3.2")
    } else { library(package, character.only = TRUE ) }
  }
}
.exactPoiCI <- function (x, conf.level = 0.95) {
  poisson.test(x, conf.level = conf.level)$conf.int[1:2]
}
.poissonCorrection <- function(posCount, count, iterations=10){
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
.mgsub <- function(pattern, replacement, x, ...) {
  if (length(pattern)!=length(replacement)) {
    stop("pattern and replacement do not have the same length.")
  }
  result <- x
  for (i in 1:length(pattern)) {
    result <- gsub(pattern[i], replacement[i], result, ...)
  }
  result
}
.Mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
getMode <- function(data, tData = NULL){
  result <- c(Mode(round(x[,1])), Mode(round(x[,2])))
  result <- thresholdData(tData = tData, amplitude = result, type = 'mode')
  return(result)
}
