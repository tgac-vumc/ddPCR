plot.ddPCR <- function(data = NULL, well = NULL, dotres = 0.7, 
                       density = 60, pch = 16, bg = "#D3D3D3", 
                       main = "ddPCR", verbose = TRUE, new = TRUE){
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
  
  if (verbose == TRUE){
    cat("Plotting sample: ", data@phenoData$sampleData['name',sample], "\n")
  }
  
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
      main.text <- paste(data@phenoData$sampleData['name',sample],":",
                         data@phenoData$sampleData['probe',sample])
      plot(0, bty='n', col = bg,
           ylab = "Ch1 Amplitude",
           xlab = "Ch2 Amplitude",
           xlim = c(0, max.ch2), 
           ylim = c(0, max.ch1),
           main = main.text
           )
      sub.text <- .dropletCountText(x = data@assayData$Cluster[,sample])
      mtext(side = 3, text = sub.text, cex = 0.8)
    
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
plot.probe <- function(data = NULL, probe = NULL, dotres = 0.7, 
                       density = 60, pch = 16, bg = "#D3D3D3", 
                       main = "ddPCR", verbose = TRUE, new = TRUE){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  } 
  probes <- unique(data@phenoData$sampleData['probe',])
  if(is.null(probe) == TRUE){
    cat("No probe was given to plot.\n")
  } 
  if(class(probe) == "integer"){
    probe <- as.numeric(probe)
  }
  if(class(probe) == "numeric"){
    if((probe > length(probes)) == TRUE){
      cat("No correct probe number has been given.\n")
    }
  } 
  if(class(density) == "numeric"){
    if(density > 0 | density < 100){
      data <- .updateColors(data = data, density = density)
    }
  }
  
  ### get highest threshold from match with sample
  selection.probe <-  data@phenoData$sampleData['probe',] %in% probes[probe]
  max.ch1 <- max(data@phenoData$ch1['maxAmplitude', selection.probe])
  max.ch2 <- max(data@phenoData$ch2['maxAmplitude', selection.probe])
  
  if (verbose == TRUE){
    cat("Plotting probe: ", probes[probe], "\n")
  }
  par(bg = bg)
  main.text <- paste("Stacked data probe:", probes[probe])
  plot(0, bty='n', col = bg,
       ylab = "Ch1 Amplitude",
       xlab = "Ch2 Amplitude",
       xlim = c(0, max.ch2), 
       ylim = c(0, max.ch1),
       main = main.text
  )
  sub.text <- .dropletCountText(x = data@assayData$Cluster[,selection.probe])
  mtext(side = 3, text = sub.text, cex = 0.8)
  
  # FOR GRID LINES
  line500 <- c(-500,0, cumsum(rep(500,50)))
  abline(v = line500, col = "#FFFFFF60")
  abline(h = line500, col = "#FFFFFF60")
  # add points
  points(y = data@assayData$Ch1.Amplitude[,selection.probe], 
         x = data@assayData$Ch2.Amplitude[,selection.probe], 
         cex = dotres, 
         col = data@assayData$Color[,selection.probe], 
         pch = pch
  )
  threshold.ch1 <- unique(data@phenoData$ch1['threshold', selection.probe], na.rm = TRUE)
  abline(h = threshold.ch1, col = "#ff0000")
  threshold.ch2 <- unique(data@phenoData$ch2['threshold', selection.probe], na.rm = TRUE)
  abline(v = threshold.ch2, col = "#ff0000")
}