Mode <- function(x){
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
getMode <- function(x, tData = NULL){
  result <- c(Mode(round(x[,1])), Mode(round(x[,2])))
  result <- thresholdData(tData = tData, amplitude = result, type = 'mode')
  return(result)
}
