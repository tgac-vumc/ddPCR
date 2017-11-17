plotThresholds <- function(thresholds, col = "black"){
  abline(h = thresholds[1], col = col)
  abline(v = thresholds[2], col = col)
}
plot.ddPCR <- function(data = NULL, dotres = 0.7, tData = NULL, main = "ddPCR", pch = 16, density = 60, verbose = FALSE){
  if(is.null(data) == TRUE){
    stop("data is not given.\n")
  }
  if(is.null(tData) != TRUE){
    if (class(tData) == "matrix") {
      if ("maxAmplitude" %in% row.names(tData) == TRUE){
        if(verbose == TRUE){cat("Reading Threshold data.\n")}
        xmax <- tData[row.names(tData) %in% "maxAmplitude",2]
        ymax <- tData[row.names(tData) %in% "maxAmplitude",1]
      } 
    } else {
      if(verbose == TRUE){cat("Amplitude data not found. tData will not be used.\n")}
      xmax <- max(data[,2])
      ymax <- max(data[,1])
    }
  } else {
    if(verbose == TRUE){cat("Threshold data must be a matrix. tData will not be used.\n")}
    xmax <- max(data[,2])
    ymax <- max(data[,1])
    }
   
  col.vec <- .defineColor(x = data[,3], density = density)
  plot(y = data[,1], x = data[,2], cex = dotres, col = col.vec, 
       ylab = "Ch1 Amplitude", xlab = "Ch2 Amplitude", 
       pch = pch, main = main,
       xlim = c(0, xmax), ylim = c(0, ymax))
  sub.text <- dropletCountText(x = data[,3])
  mtext(side = 3, text = sub.text, cex = 0.8)
  
  if (class(tData) != "matrix") {
    if(verbose == TRUE){cat("Threshold data must be a matrix. Data will not be plotted.")}
  } else {
    if("thresholdMeanStDev" %in% row.names(tData) == TRUE) 
      {
      abline(v = tData[row.names(tData) %in% "thresholdMeanStDev",2], col = "red")
      abline(h = tData[row.names(tData) %in% "thresholdMeanStDev",1], col = "red")
      }else if("threshold" %in% row.names(tData) == TRUE) 
        {
      abline(v = tData[row.names(tData) %in% "threshold",2], col = "red")
      abline(h = tData[row.names(tData) %in% "threshold",1], col = "red")
        }
    }
}