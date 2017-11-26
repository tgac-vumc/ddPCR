ddPCRdata <- setClass("ddPCRdata", slots = c(assayData = "list", 
                                             phenoData = "list", 
                                             experimentData = "list"))
readAmplitudeFiles <- function(files =  NULL, nrows=25000, verbose=FALSE){
  if(is.null(files) == TRUE){
    stop("No files input has been given.\n\n")
  } else if(FALSE %in% file.exists(files) == TRUE){
    stop("Not all files are present.\n\n")
  } else {
    sample.names <- .mgsub(x=basename(files), pattern = c("_Amplitude.csv",".csv"), replacement = c("",""))
    well.locations <- .findWell(basename(files))
  }
  nFiles <- length(files)
  
  ### CREATE BASE MATRIX FOR ALL DATA ------
  Ch1.Amplitude <- matrix(NA, nrow = nrows, ncol = nFiles, dimnames = list(c(1:nrows), c(well.locations)))
  Ch2.Amplitude <- matrix(NA, nrow = nrows, ncol = nFiles, dimnames = list(c(1:nrows), c(well.locations)))
  Cluster <- matrix(NA, nrow = nrows, ncol = nFiles, dimnames = list(c(1:nrows), c(well.locations)))
  Color <- matrix(NA, nrow = nrows, ncol = nFiles, dimnames = list(c(1:nrows), c(well.locations)))
  Use <- matrix(FALSE, nrow = nrows, ncol = nFiles, dimnames = list(c(1:nrows), c(well.locations)))
  
  ### ADD PHENODATA  -----------
  rowsSampledata <- c("name", "file", "well", "probe", "type")
  sample.data <- matrix(data = NA, nrow = length(rowsSampledata), 
                        ncol =  length(sample.names),
                        dimnames = list(rowsSampledata, well.locations))
    
  rowsPhenodata <- c("totalDroplets","minOutlier","maxOutlier","threshold","minAmplitude","maxAmplitude","minRain","maxRain")
  channel.1 <- matrix(data = NA, nrow = length(rowsPhenodata), 
                      ncol =  length(sample.names),
                      dimnames = list(rowsPhenodata, well.locations))
  channel.2 <- matrix(data = NA, nrow = length(rowsPhenodata), 
                      ncol =  length(sample.names),
                      dimnames = list(rowsPhenodata, well.locations))
  
  Channel.Name <- c("Ch1.Amplitude", "Ch2.Amplitude") 
  Cluster.Number <- c(0, 1, 2, 3, 4, 5)
  Cluster.Name <- c("outlier", "ch1-ch2-", "ch1+ch2-", "ch1+ch2+", "ch1-ch+", "rain")
  Cluster.Color <- c("#FF0000", "#000000", "#0000FF", "#FFA500", "#00FF00", "#551A8B")
  
  ### loop through all the files ------
  for(i in 1:nFiles)
  {
    if(file.exists(files[i]) == TRUE)
    {
      amplitude.data <- read.table(file = files[i], header = TRUE,sep = ",")
      if(verbose == TRUE){ 
        if(dim(amplitude.data)[1] > nrows){
          cat("Too many rows in data file: '", basename(files[i]),"'\n", sep="")
          next
        } else {cat("Reading amplitude data: '", basename(files[i]),"'\n", sep="")}
      }
      
      Total.Droplets <- nrow(amplitude.data)
      ### Add sample data -----------
      Ch1.Amplitude[1:Total.Droplets,i] <- amplitude.data[,1]
      Ch2.Amplitude[1:Total.Droplets,i] <- amplitude.data[,2]
      Cluster[1:Total.Droplets,i] <- amplitude.data[,3]
      Use[1:Total.Droplets,i] <- TRUE
      
      channel.1['totalDroplets', i] <- nrow(amplitude.data) # TotalDroplets
      channel.2['totalDroplets', i] <- nrow(amplitude.data) # TotalDroplets
      
      channel.1['minAmplitude', i] <- min(amplitude.data[,1])
      channel.1['maxAmplitude', i] <- max(amplitude.data[,1])
      channel.2['minAmplitude', i] <- min(amplitude.data[,2])
      channel.2['maxAmplitude', i] <- max(amplitude.data[,2])
    }
  }
  
  ### ADD SAMPLE.DATA TO STRUCTURE ------
  sample.data['well',] <- well.locations
  sample.data['name',] <- sample.names
  sample.data['file', ] <- basename(files)
  
  ### read experiment 
  design <- .readExperiment(files, verbose = verbose)
  matching <- match(sample.data['file',], design[,2])
  sample.data['name',] <- design[matching,1]
  sample.data['probe',] <- design[matching,4]
  sample.data['type',] <- design[matching,3]
  

  
  ### COMBINE ALL DATA INTO ONE STRUCTURE ------
  assayData <- list(Ch1.Amplitude = Ch1.Amplitude, 
                    Ch2.Amplitude = Ch2.Amplitude, 
                    Cluster = Cluster, 
                    Use = Use, 
                    Color = Color)
  phenoData <- list(sampleData = sample.data,
                    ch1 = channel.1, 
                    ch2 = channel.2)
  experimentData <- list(Cluster.Number = Cluster.Number,
                         Cluster.Name = Cluster.Name,
                         Cluster.Color = Cluster.Color
  )
  
  object <- new('ddPCRdata', assayData=assayData, phenoData=phenoData, 
                experimentData=experimentData)
  object <- .updateColors(data = object)
  
  return(object)
}