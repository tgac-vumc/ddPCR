plot_ddPCR <- function(x, dotres = 0.7, main = "ddPCR", pch = 16, colors = "ddpcr", density = 60, thresholds = NULL, max.xy = NULL, verbose = FALSE)
{
  if(length(max.xy) != 2) {
    xmax <- max(x[,2])
    ymax <- max(x[,1])
  } else {
    xmax <- max.xy[2]
    ymax <- max.xy[1]
  }
  col.vec <- defineColor(x = x[,3], density = density)
  plot(y = x[,1], x = x[,2], cex = dotres, col = col.vec, 
       ylab = "Ch1 Amplitude", xlab = "Ch2 Amplitude", 
       pch = pch, main = main,
       xlim = c(0, xmax), ylim = c(0, ymax))
  sub.text <- dropletCountText(x = x[,3])
  mtext(side = 3,text = sub.text, cex = 0.8)
  if (class(thresholds) != "matrix") {
    if(verbose == TRUE){cat("Threshold data must be a matrix. Data will not be plotted.")}
  }else{
    abline(h = thresholds[,1], col = "red") # channel 1
    abline(v = thresholds[,2], col = "red") # channel 2
  }
}