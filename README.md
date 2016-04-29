# README file for ddPCR
# 
# 201604, HFvanEssen

Use ddPCR QuantaSoft software to export raw Amplitude data
- 'Export Amplitude and Cluster Data' in options 

Pipeline example:
- create design
- loop through individual targets
	- combine all samples for target 
	- find global thresholds combined data
	- refine global thresholds for combined data
	- define clusters for combined data
	- plot combined data with refined thresholds
	- loop through each sample for each target with refined thresholds
	- plot data
	- use calculateStats... functions to get sample statistics


######
anaysis.path <- "set your path here..."

createDesign(path = analysis.path)
design <- read.table(file.path(path,"design.txt"), header = TRUE, sep = ",")

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
      data.xy.max <-  getMaxAmplitude(combined.data)
      thresholds <- getThresholds(combined.data, algorithm = "hist")
      combined.data <- defineClusters(combined.data, thresholds)
      
      sample.name <- "Combined_data_base_threshold"
      output.file <- file.path(path.targets[i], paste(sample.name,".png",sep=""))
      png(filename = output.file, width = 800, height = 800)
      plot_ddpcr(x = combined.data, main = sample.name, max.xy = data.xy.max, thresholds = thresholds)
      dev.off()
      
      thresholds.2 <- 
        combined.data %>% 
        refineThresholdStDev(., stdev = 3, thresholds = thresholds)
      data.threshold.2 <- defineClusters(combined.data, thresholds.2)
      
      sample.name <- "Combined_data_base_threshold_2"
      output.file <- file.path(path.targets[i], paste(sample.name,".png",sep=""))
      png(filename=output.file, width = 800, height = 800)
      plot_ddpcr(x=data.threshold.2, main=sample.name, max.xy=data.xy.max, thresholds=thresholds.2)
      dev.off()
   
      results <- c()
      for(j in 1:length(files))
      {
        sample.data <- 
          read.table(file=file.path(path, files[j]),header = TRUE,sep = ",") %>%
          defineClusters(., thresholds.2)
        output.file <- file.path(path.targets[i], paste(file.wells[j],"_",file.names[j],".png",sep=""))
        png(filename=output.file,width = 800,height = 800)
        plot_ddpcr(x=sample.data, main=paste(file.wells[j],": ",file.names[j]), max.xy = data.xy.max, thresholds = thresholds.2)
        dev.off()
        
        result <- NULL
        result <- cbind(Well=rep(as.character(file.wells[j]),2),Sample=rep(as.character(file.names[j]),2))
        result <- data.frame(result)
        result <- cbind(result, TargetType=c("Channel 1","Channel 2"))
        result <- cbind(result, Target=rep(targets[i],2))
        result <- cbind(result, Status=rep(sampleQC(x = sample.data, sample.type = sample.type[j]),2))
        result <- cbind(result, Threshold=thresholds.2)
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

ddpcr.analysis(path = analysis.path, experimentDesign=design)

