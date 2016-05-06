library(tcltk)

### DDPCR ANALYSIS PIPELINE
ddpcr.analysis <- function(path, experimentDesign)
{
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
    
    if(sum(files %in% list.files(path)) == length(files))
    {
      combined.data <- combineSamples(path = path, files = files) 
      tData <- defineMinOutliers(combined.data)
      tData <- getMaxAmplitude(combined.data, tData = tData)
      tData <- getThresholds(combined.data, algorithm = "hist", tData = tData)
      combined.data <- defineClusters(combined.data, thresholds)
      
      sample.name <- "Combined_data_base_threshold"
      output.file <- file.path(path.targets[i], paste(sample.name,".png",sep=""))
      png(filename = output.file, width = 800, height = 800)
      plot_ddPCR(x = combined.data, main = sample.name, max.xy = data.xy.max, thresholds = thresholds)
      dev.off()
      
      thresholds.2 <- refineThresholdStDev(combined.data, stdev = 3, thresholds = thresholds)
      data.threshold.2 <- defineClusters(combined.data, thresholds.2)
      
      sample.name <- "Combined_data_base_threshold_2"
      output.file <- file.path(path.targets[i], paste(sample.name,".png",sep=""))
      png(filename=output.file, width = 800, height = 800)
      plot_ddPCR(x=data.threshold.2, main=sample.name, max.xy=data.xy.max, thresholds=thresholds.2)
      dev.off()
   
      results <- c()
      for(j in 1:length(files))
      {
        sample.data <- 
          read.table(file=file.path(path, files[j]),header = TRUE,sep = ",") %>%
          defineClusters(., thresholds.2)
        output.file <- file.path(path.targets[i], paste(file.wells[j],"_",file.names[j],".png",sep=""))
        png(filename=output.file,width = 800,height = 800)
        plot_ddPCR(x=sample.data, main=paste(file.wells[j],": ",file.names[j]), max.xy = data.xy.max, thresholds = thresholds.2)
        dev.off()
        
        result <- NULL
        result <- cbind(Well=rep(as.character(file.wells[j]),2),Sample=rep(as.character(file.names[j]),2))
        result <- data.frame(result)
        result <- cbind(result, TargetType=c("Channel 1","Channel 2"))
        result <- cbind(result, Target=rep(targets[i],2))
        result <- cbind(result, Status=rep(sampleQC(x = sample.data, sample.type = sample.type[j]),2))
        # result <- cbind(result, Threshold=thresholds.2)
        result <- cbind(result, calculateStatsDroplets(sample.data))
        copies.data <- CalculateStatsCopies(sample.data)
        result <- cbind(result, copies.data)
        result <- cbind(result, ngPer1ul=convertCopiesToNanogram(result$CopiesPer1ul))
        result <- cbind(result, calculateStatsRatioFract(copies.data))
        results <- rbind(results,result)
      }
      
      output.file <- file.path(path.targets[i],paste(experiment,"_",targets[i], "_results.txt",sep=""))
      write.table(x = results,file = output.file,quote = FALSE,sep = "\t",row.names = FALSE)
      cat("Probe", as.character(targets[i]), "has been processed.\n")
    } else {cat("Not all files for target ",targets[[i]]," are present for processing.\n", sep="")}
  }
}

####
# SCRIPT USES A DESIGN FILE TO RUN THE PIPELINE
####


### SOURCE FUNCTIONS AND SET ANALYSIS PATH
computer.name <- as.character(Sys.info()["nodename"]) 
if(computer.name == "localhost.localdomain")
{
  source.path <- "/home/dirk/Documents/r_scripts/ddpcr_analysis/scripts_ddpcr"
  R.utils::sourceDirectory(source.path, modifiedOnly=FALSE);
  analysis.path <- "/run/user/1000/gvfs/smb-share:server=vumc-cl-fs02-01,share=microarray-faciliteit$/Data ddPCR/20160418 EGFR spike-in test sm_2016-04-18-15-19"
}
if(computer.name == "MacBook-Air-van-Dirk.local")
{ 
  source.path <- ""
  R.utils::sourceDirectory(source.path, modifiedOnly=FALSE);
  analysis.path <- "/Users/dirkvanessen/ownCloud/r_scripts/ddPCR analysis/input.data/20160317 EGFR test sm_2016-03-17-14-59"
}
### CREATE DESIGN
if(file.exists(file.path(analysis.path, "design.txt")) == FALSE)
{
  createDesign(path=analysis.path)
}
### LOAD DESIGN / RUN ANALYSIS
test <- .rs.askForPassword(prompt = "Is design file ready (y/n)")
if(tolower(test) %in% c("y","yes","true") == TRUE)
{
  experimentDesign <- read.table(file.path(analysis.path,"design.txt"), header = TRUE, sep = ",")
  cat ("Processing ddPCR run.\n")
  ddpcr.analysis(path = analysis.path, experimentDesign = experimentDesign)
  
}


  


