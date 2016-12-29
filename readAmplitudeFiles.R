ddPCRdata <- setClass("ddPCRdata", slots = c(assayData = "list", phenoData = "matrix", experimentData = "list"))

readAmplitudeFiles <- function(files, nrows=25000, verbose=FALSE){ # not used in pipeline
  nFiles <- length(files)
  Ch1.Amplitude <- matrix(NA, nrow = nrows, ncol = nFiles)
  Ch2.Amplitude <- matrix(NA, nrow = nrows, ncol = nFiles)
  Cluster <- matrix(NA, nrow = nrows, ncol = nFiles)
  Color <- matrix(NA, nrow = nrows, ncol = nFiles)
  Use <- matrix(FALSE, nrow = nrows, ncol = nFiles)
  Total.Droplets <- rep(NA, nFiles)
  Sample.Names <- rep("NA", nFiles)
  Ch1.Max <- rep(NA, nFiles)
  Ch2.Max <- rep(NA, nFiles)
  
  for(i in 1:nFiles)
  {
    if(file.exists(files[i]) == TRUE)
    {
      sample.data <- read.table(file = files[i], header = TRUE,sep = ",")
      if(verbose == TRUE){ 
        if(dim(sample.data)[1] > nrows){
          cat("Too many rows in data file: '", basename(files[i]),"'\n", sep="")
          next
        } else {cat("Reading amplitude data: '", basename(files[i]),"'\n", sep="")}
      }
      Sample.Names[i] <- mgsub(x=basename(files[i]), pattern = c("_Amplitude.csv",".csv"), replacement = c("",""))
      Total.Droplets[i] <- dim(sample.data)[1]
      Ch1.Amplitude[1:Total.Droplets[i],i] <- sample.data[,1]
      Ch2.Amplitude[1:Total.Droplets[i],i] <- sample.data[,2]
      Cluster[1:Total.Droplets[i],i] <- sample.data[,3]
      Use[1:Total.Droplets[i],i] <- TRUE
      Ch1.Max[i] <- max(sample.data[,1])
      Ch2.Max[i] <- max(sample.data[,2])
    }
  }
  assayData <- list(Ch1.Amplitude = Ch1.Amplitude, Ch2.Amplitude = Ch2.Amplitude, Cluster = Cluster, Use = Use)
  phenoData <- cbind(names=as.character(Sample.Names), Total.Droplets=Total.Droplets, Ch1.Max=Ch1.Max, Ch2.Max=Ch2.Max)
  experimentData <- list(
    Channel.Name=c("Ch1.Amplitude" ,"Ch2.Amplitude"), 
    Cluster.Number=c(1,2,3,4,0,5),
    Cluster.Name=c("ch1-ch2-", "ch1+ch2-", "ch1+ch2+", "ch1-ch+", "outlier", "rain"),
    Cluster.Color=c("black", "blue", "orange", "green", "red", "purple")
  )
  object <- new('ddPCRdata', assayData=assayData, phenoData=phenoData, experimentData=experimentData)
  return(object)
}
combineSamples <- function(path, files, verbose = FALSE){
  combined.data <- NULL
  if(verbose == TRUE){
    cat("Reading Amplitude files.\n")
  }
  for(i in 1:length(files))
  { 
    nrlines <- NULL
    file.name <- file.path(path, files[i])
      if(length(count.fields(file.name)) > 0){
      sample.data <- read.table(file=file.name, header = TRUE, sep = ",")
      combined.data <- rbind(combined.data, sample.data )
      }
  }
  if(verbose == TRUE){
    cat("Amplitude files are combined.\n")
  }
  return(combined.data)
}