plot.ddPCR <- function(data = NULL, well = NULL, dotres = 0.7, density = 60, pch = 16, bg = "lightgrey", main = "ddPCR", verbose = FALSE, test = FALSE){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  } 
  if(is.null(well) == TRUE){
    cat("No well location was given to plot.\n")
  } else {
    sample <- match(tolower(well), table = tolower(data@experimentData$Sample.Well))
    if(is.na(sample) == TRUE){
      cat("No correct sample well has been given.\n")
    } else if(test == FALSE){
      
      data <- .updateColors(data = data, density = density)
      
      plot(y = data@assayData$Ch1.Amplitude[,sample], 
           x = data@assayData$Ch2.Amplitude[,sample], 
           cex = dotres, 
           col = data@assayData$Color[,sample], 
           xlab = "Ch2 Amplitude", 
           ylab = "Ch1 Amplitude",
           pch = pch, 
           main = toupper(well),
           xlim = c(0, max(data@phenoData$ch2[, "maxAmplitude"])), 
           ylim = c(0, max(data@phenoData$ch1[, "maxAmplitude"])))
      sub.text <- .dropletCountText(x = data@assayData$Cluster[,sample])
      mtext(side = 3, text = sub.text, cex = 0.8)
    } else {
      par(bg = bg)
      plot(0, bty='n', col = bg,
           ylab = "Ch1 Amplitude",
           xlab = "Ch2 Amplitude",
           xlim = c(0, max(data@phenoData$ch2[, "maxAmplitude"])), 
           ylim = c(0, max(data@phenoData$ch1[, "maxAmplitude"])),
           main = toupper(well)
           )
    
      # FOR GRID LINES
      line500 <- c(-500,0, cumsum(rep(500,50)))
      abline(v = line500, col = "#FFFFFF60")
      abline(h = line500, col = "#FFFFFF60")
      # add points
      points(y = data@assayData$Ch1.Amplitude[,sample], 
             x = data@assayData$Ch2.Amplitude[,sample], 
             cex = dotres, 
             col = data@assayData$Color[,sample], 
             pch = pch
      )
    }
  }
}