createDesign <- function(path, verbose = TRUE){ 
  if(!file.exists(file.path(path, "design.txt"))) {
  experiment <- list.files(path, pattern = "Amplitude.csv", full.names = FALSE)[1]
  experiment <- strsplit(x = experiment, split = "_")[[1]][1]
  experiment.file <- paste(experiment, ".csv", sep = "")

  if(file.exists(file.path(path, experiment.file)) == TRUE){
    if(verbose == TRUE){
      cat("Experiment setup has been located.\n")
      }
    data <- read.table(file = file.path(path, experiment.file), header = TRUE, sep = ",", check.names = FALSE, row.names = NULL)
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
      if(tolower(sampleName) %in% c("ntc", "water", "neg", "negative", "te buffer", "h2o") == TRUE)
      {
        design[well, 3] <- "neg"
      } else { design[well, 3] <- "sample" }
      design[well, 4] <- paste(as.character(sort(temp$Target)), collapse = " vs ")
    }
  } else  
  {
    if(verbose == TRUE){
      cat("Experiment setup has not been located. Basic experiment design file will be made.\n")
      }
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
getControls <- function(x, pos = c("positive", "pos control"), ntc = c("te buffer", "water"), neg = ""){ 
  # x = vector of sample names is input
  results <- rep("sample", length(x))
  x <- tolower(x)
  results[x %in% tolower(pos)] <- "pos"
  results[x %in% tolower(ntc)] <- "ntc"
  results[x %in% tolower(neg)] <- "neg"
  return(results)
}
getTargets <- function(path){ # will use .csv file of the experiment located in path
  experiment.file <- list.files(path, pattern = "Amplitude.csv", full.names = FALSE)[1]
  experiment.file <- paste(strsplit(x = experiment.file, split = "_")[[1]][1], ".csv", sep = "")
  if(file.exists(experiment.file) == TRUE)
  {
    data <-  read.csv(file = file, header = TRUE, sep = ",", check.names = FALSE, row.names = NULL)
    if(colnames(data)[1] == "row.names")
    {
      colnames(data) <- colnames(data)[2:dim(data)[2]]
      data <- data[,1:(dim(data)[2]-1)]
    }
  } else { break 
  }
  targets <- unique(mgsub(pattern = c("_WT", "_wt", " wt", " WT"), replacement = c("", "", "", ""), x = data$Target))
  result <- list()
  for(i in 1:length(targets))
  {
    result[[i]] <- data[grep(data$Target,pattern = targets[i]),]
  }
  targets <- gsub(pattern = " ", replacement = "_", x = targets)
  names(result) <- targets
  return(result)
}
createTargetFolders <- function(x, path = NULL){ # x = output of get.targets(). Will create output folders at specified path.
  paths <- NULL
  if(class(x) == "list" & class(path) != "NULL")
  {
    paths <- file.path(path, names(x))
    for(i in 1:length(paths))
    {
      if(!file.exists(paths[i]))
      {
        dir.create(paths[i])
      }
    }
  } else if (class(x) == "character" | class(x) == "factor")
  {
    if(class(path) != "NULL")
    {
      for(i in 1:length(x))
      {
        paths[i] <- file.path(path, x[i])
        if(!file.exists(paths[i]))
        {
          dir.create(paths[i])
        }
      }
    }
  }
  return(paths)
}
# [ ] getTargets has overlap with createDesign... 
# [ ] getControls has overlap with createDesign... 