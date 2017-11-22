.getExperiment <- function(file){
  file <- basename(file)
  file.length <- nchar(file)
  experiment <- substr(x = file,  start = 1, stop = (file.length-18))
  return(experiment)
}
.getPath <- function(file = NULL){
  if(is.null(file) == TRUE){
    stop("No file has been given.\n")
  } else if(length(files) > 1){
      file <- file[1]
    }
  basenameFile <- basename(file)
  path <- gsub(file, pattern = basenameFile, replacement = "")
  
  return(path)
}
.readExperiment <- function(file = NULL, verbose = TRUE){
  path <- .getPath(file)
  experiment <- unique(getExperiment(basename(file)))
  experiment.file <- paste0(experiment, ".csv")
  if(file.exists(file.path(path,experiment.file)) == TRUE){
    if(verbose == TRUE){
      if(length(experiment) > 1){
      stop("Multiple Experiment setups have been located.\n Not working with this type of setup yet.")
    } else { 
      cat("Experiment setup has been located.\n")
     }
    }
  }
  experimentData <- NULL
  #for(i in 1:length(experiment)){
    data <- read.table(file = file.path(path, experiment.file), header = TRUE, sep = ",", check.names = FALSE, row.names = NULL, fill = TRUE)
    experimentData <- rbind(experimentData, data)
    if(ncol(data) == 1){
      if(verbose == TRUE){
        cat("File is not separated by ',' \nChecking alternative.....")
      }
      data <- read.table(file = file.path(path, experiment.file.file), header = TRUE, sep = ";", check.names = FALSE, row.names = NULL, fill = TRUE)
      experimentData <- rbind(experimentData, data)
      if(ncol(data) == 1){
        stop("Data is not provided in the correct format.\n")
      }
      if(verbose == TRUE){
        cat("     Alternative accepted.\n")
      }
    }
  #}
    #### create design
    data <- data.frame(experimentData[,1:9])
    colnames(data) <- c("Well","AmplitudeFile","Experiment","Sample","TargetType","Target","Status","Concentration","Supermix")
    wells <- unique(data$Well)
    data$AmplitudeFile <-  paste(experiment, "_", data$Well, "_Amplitude.csv", sep = "")
    amplitude.files <- list.files(path = path, pattern = "_Amplitude.csv")
    data$AmplitudeFile <- amplitude.files[match(x = data$AmplitudeFile, table = amplitude.files)]
    data <- data[!is.na(data$AmplitudeFile),]
 
    design <- matrix(data = "", nrow = length(wells), ncol = 4)
    colnames(design) <- c("Name", "File", "Type", "Probe")
    for(well in 1:length(wells)) {
      temp <- data[data$Well %in% wells[well],]
      sampleName <- unique(temp$Sample)
      design[well, 1] <- paste(sampleName, collapse = "_")
      design[well, 2] <- paste(experiment, "_", wells[well], "_Amplitude.csv", sep = "")
      if(tolower(sampleName) %in% c("ntc", "water", "neg", "negative", "te buffer", "h2o") == TRUE)
      {
        design[well, 3] <- "neg"
      } else { design[well, 3] <- "sample" }
      design[well, 4] <- paste(as.character(sort(temp$Target)), collapse = " vs ")
    }
    
    return(design)
}

