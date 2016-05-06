createDesign <- function(path, verbose = TRUE)
{ if(!file.exists(file.path(path, "design.txt")))
  {
  experiment <- list.files(path, pattern = "Error.log", full.names = FALSE)
  experiment <- gsub(pattern = "Error.log", replacement = "", x = experiment)
  experiment.file <- file.path(path, paste(experiment, ".csv", sep = ""))
  if(file.exists(experiment.file) == TRUE)
  {
    if(verbose == TRUE){cat("Experiment setup has been located.\n")}
    data <- read.table(file = experiment.file, header = TRUE, sep = ",", check.names = FALSE, row.names = NULL)
    colnames <- c(colnames(data)[2:9],"Supermix")
    data <- data.frame(data[,1:9])
    colnames(data) <- colnames
    data$ExptType <-  paste(experiment, "_", data$Well, "_Amplitude.csv", sep = "")
    colnames(data)[2] <- "AmplitudeFile"
    amplitude.files <- list.files(path = path, pattern = "_Amplitude.csv")
    data <- data[data$AmplitudeFile %in% amplitude.files,]
    
    wells <- unique(data$Well)
    
    design <- matrix(data = "", nrow = length(wells), ncol = 4)
    colnames(design) <- c("Name", "File", "Type", "Probe")
    
    for(well in 1:length(wells))
    {
      temp <- data[data$Well %in% wells[well],]
      sampleName <- unique(temp$Sample)
      design[well, 1] <- paste(sampleName, collapse = "_")
      design[well, 2] <- paste(experiment, "_", wells[well], "_Amplitude.csv", sep = "")
      if(tolower(sampleName) %in% c("ntc", "water", "neg", "negative") == TRUE)
      {
        design[well, 3] <- "neg"
      } else { design[well, 3] <- "sample" }
      design[well, 4] <- paste(as.character(sort(temp$Target)), collapse = " vs ")
    }
  } else  
  {
    if(verbose == TRUE){cat("Experiment setup has not been located. Basic experiment design file will be made.\n")}
    amplitude.files <- list.files(path = path, pattern = "_Amplitude.csv")
    sample.names <- gsub(pattern = "_Amplitude.csv", x = amplitude.files, replacement = "")
    design <- matrix(data = "", nrow = length(amplitude.files), ncol = 4)
    colnames(design) <- c("Name", "File", "Type", "Probe")
    design[,1] <- sample.names
    design[,2] <- amplitude.files
    design[,3] <- c("pos", (rep("sample", length(amplitude.files)-2)), "neg")
    design[,4] <- "probe"
  }
  output.file <- file.path(path, "design.txt")
  write.table(file = output.file, x = design, quote = FALSE, sep = "\t", row.names = FALSE)
  }
}
