plot.ddPCR <- function(data = NULL, well = NULL, dotres = 0.7, 
                       density = 60, pch = 16, bg = "#e6e6e6", 
                       main = "ddPCR", verbose = FALSE, new = FALSE){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  } 

  if(is.null(well) == TRUE){
    cat("No well location was given to plot.\n")
  } 
  if(class(well) == "integer"){
    well <- as.numeric(well)
  }
  if (class(well) == "character"){
    sample <- match(tolower(well), table = tolower(data@phenoData$sampleData['well',]))
    if(is.na(sample) == TRUE){
        cat("No correct sample well has been given.\n")
    }
  } else if(class(well) == "numeric"){
      if((well > ncol(data@phenoData$sampleData)) == TRUE){
        cat("No correct sample well has been given.\n")
      }else {
        sample <- well
      }
  } 
  if(class(density) == "numeric"){
    if(density > 0 | density < 100){
      data <- .updateColors(data = data, density = density)
    }
  }

  ### get highest threshold from match with sample
  selection.probe <-  data@phenoData$sampleData['probe',] %in% data@phenoData$sampleData['probe',sample]
  max.ch1 <- max(data@phenoData$ch1['maxAmplitude', selection.probe])
  max.ch2 <- max(data@phenoData$ch2['maxAmplitude', selection.probe])
  
  cat("working on sample: ", data@phenoData$sampleData['name',sample], "\n")
   if(new != TRUE){
      data <- .updateColors(data = data, density = density)
    
      plot(y = data@assayData$Ch1.Amplitude[,sample], 
           x = data@assayData$Ch2.Amplitude[,sample], 
           cex = dotres, 
           col = data@assayData$Color[,sample], 
           xlab = "Ch2 Amplitude", 
           ylab = "Ch1 Amplitude",
           pch = pch, 
           main = data@phenoData$sampleData['name', sample],
           xlim = c(0, max.ch2), 
           ylim = c(0, max.ch1))
      sub.text <- .dropletCountText(x = data@assayData$Cluster[,sample])
      mtext(side = 3, text = sub.text, cex = 0.8)
      
      if(is.na(data@phenoData$ch1['threshold', sample]) != TRUE){
        abline(h = data@phenoData$ch1['threshold', sample], col = "#ff0000")
      }
      if(is.na(data@phenoData$ch2['threshold', sample]) != TRUE){
        abline(v = data@phenoData$ch2['threshold', sample], col = "#ff0000")
      }
    } else if (new == TRUE) {
      par(bg = bg)
      plot(0, bty='n', col = bg,
           ylab = "Ch1 Amplitude",
           xlab = "Ch2 Amplitude",
           xlim = c(0, max.ch2), 
           ylim = c(0, max.ch1),
           main = data@phenoData$sampleData['name',sample]
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
      if(is.na(data@phenoData$ch1['threshold', sample]) != TRUE){
        abline(h = data@phenoData$ch1['threshold', sample], col = "#ff0000")
      }
      if(is.na(data@phenoData$ch2['threshold', sample]) != TRUE){
        abline(v = data@phenoData$ch2['threshold', sample], col = "#ff0000")
      }
    }
  }
