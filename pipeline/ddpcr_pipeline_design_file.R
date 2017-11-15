# This is and example of the ddPCR pipeline analysis.
# Before use, the file MUST be copied to a different
# location. Then the "#" may be removed and the code
# can be run. Please do not make any changes in this
# document.
#
# HFvE, 20161227

#### START USER INPUT

source.path <- "enter path with location of ddPCR functions here..."

#### END USER INPUT

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,"")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
ddpcr.analysis <- function(path) {
  createDesign(path = path)
  experimentDesign <- read.table(file.path(path,"design.txt"), header = TRUE, sep = "\t")
  targets <- unique(experimentDesign$Probe)
  path.targets <- createTargetFolders(x = targets, path = path)
  
  for(i in 1:length(targets))
  {
    experiment <- basename(path)
    design <- experimentDesign[experimentDesign$Probe == targets[i],]
    files  <- design$File
    sample.type <- design$Type
    file.names <- design$Name
    file.wells <- sapply(X = design$File, FUN = findWell)
    
    if(sum(files %in% list.files(path)) == length(files)) # combined data
    {
      # - [x] get max Amplitude of all the files
      combined.data <- combineSamples(path = path, files = files)
      tData <- getMaxAmplitude(combined.data)
      tData <- defineMinOutliers(combined.data, tData = tData)
      tData <- getThresholds(combined.data, algorithm = "hist", rm.outliers = TRUE, tData = tData)
      combined.data <- defineClusters(combined.data, tData = tData)
      
      fileName <- "combined_data_threshold"
      output.file <- file.path(path.targets[i], paste(fileName,".png",sep=""))
      png(filename=output.file,width = 800, height = 800)
      plot_ddPCR(combined.data, tData = tData, main=fileName)
      dev.off()
      
      tData <- refineThresholdStDev(combined.data, tData = tData)
      combined.data.2 <- defineClusters(combined.data, tData)
      fileName <- "combined_data_refinedThreshold"
      output.file <- file.path(path.targets[i], paste(fileName,".png",sep=""))
      png(filename=output.file,width = 800, height = 800)
      plot_ddPCR(combined.data.2, tData = tData, main = fileName)
      dev.off()
    }
      results <- NULL
      for(j in 1:length(files)) # each sample separately
      {
        sample.data <- read.table(file=file.path(path,files[j]), header = TRUE, sep = ",")
        sample.data <- defineClusters(x = sample.data, tData = tData)
        
        fileName <- paste(file.wells[j],"_", file.names[j], sep = "")
        output.file <- file.path(path.targets[i], paste(fileName, ".png", sep=""))
        png(filename = output.file, width = 800, height = 800)
        plot_ddPCR(sample.data, tData = tData, main = fileName)
        dev.off()
      
        result <- cbind(Well = rep(as.character(file.wells[j]), 2), Sample = rep(as.character(file.names[j]), 2))
        result <- data.frame(result)
        result <- cbind(result, TargetType = c("Channel 1", "Channel 2"))
        result <- cbind(result, Target = rep(targets[i],2))
        result <- cbind(result, Status = rep(sampleQC(x = sample.data, sample.type = sample.type[j]), 2))
        # result <- cbind(result, Threshold=thresholds.2) # still need to change
        result <- cbind(result, calculateStatsDroplets(sample.data))
        copies.data <- CalculateStatsCopies(sample.data)
        result <- cbind(result, copies.data)
        result <- cbind(result, ngPer1ul = convertCopiesToNanogram(result$CopiesPer1ul))
        result <- cbind(result, calculateStatsRatioFract(copies.data))
        results <- rbind(results,result)
      }
      output.file <- file.path(path.targets[i], paste(experiment,"_", targets[i], "_results.txt", sep=""))
      write.table(x = results,file = output.file, quote = FALSE, sep = "\t", row.names = FALSE)
      library(xlsx)
      output.file <- file.path(path.targets[i], paste(experiment,"_", targets[i], "_results.xlsx", sep=""))
      write.xlsx(x = results, ExcelFileName = output.file, row.names = FALSE)
      cat("Probe", as.character(targets[i]), "has been processed.\n")
    } 
}

### TO RUN.. USE CODE BELOW
#
# sourceDir(source.path) 
#
# location <- file.choose()
# inputFile <- basename(location)
# pad <- unlist(strsplit(location, inputFile))
# setwd(pad)
#
# ddpcr.analysis(path = pad)

