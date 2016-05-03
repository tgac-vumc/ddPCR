library(tcltk)
ddpcr.analysis <- function(path)
{
  experimentDesign <- read.table(file.path(path,"design.txt"), header = TRUE, sep = "\t")
  # experimentDesign  <- read.table(file=file.path(data.path,"design.txt"), header=TRUE, sep = "\t")
  targets <- unique(experimentDesign$Probe)
  # - [x] create new directory for the target
  path.targets <- createTargetFolders(x = targets, path = path)
  
  # - [x] loop through all the targets
  for(i in 1:length(targets))
  {
    experiment <- basename(path)
    design <- experimentDesign[experimentDesign$Probe == targets[i],]
    files  <- design$File
    sample.type <- design$Type
    file.names <- design$Name
    file.wells <- sapply(X = design$File, FUN = findWell)
    
    #  sample.list  <- data.targets[[i]])
    
    if(sum(files %in% list.files(path)) == length(files))
    {# are all files available for analysis.
      # - [x] get max Amplitude of all the files
      data.xy.max <- 
        combineSamples(path=path,files=files) %>%
        getMaxAmplitude(.)
      # - [x] get positive sample and determine threshold
      control.data.pos  <- 
        files[sample.type == "pos"] %>%
        combineSamples(path=path,files=.)
      # - [x] get positive control threshold
      thresholds <- 
        control.data.pos %>% 
        getThresholds(., algorithm = "hist")
      # - [x] set clusters positive control with thresholds
      control.data.pos %<>%
        defineClusters(., thresholds)
      # - [x] set file name control sample 
      control.name <- paste(file.names[sample.type == "pos"],"_pos_Control",sep="", collapse="")
      output.file <- file.path(path.targets[i], paste(control.name,".png",sep=""))
      # - [x] create plot for control data
      png(filename=output.file,width = 800,height = 800)
      plot.ddpcr(x=control.data.pos, main=control.name, max.xy=data.xy.max, thresholds=thresholds)
      dev.off()
      
      # - [x] get ntc sample(s) 
      control.data.ntc <- 
        files[sample.type == "ntc" | sample.type == "neg"] %>%
        combineSamples(path=path,files=.)
      # - [ ] set clusters ntc control with thresholds
      control.data.ntc %<>%
        defineClusters(., thresholds)
      # - [x] redefine clusters with mean & stdev
      thresholds.2 <- 
        control.data.ntc %>% 
        refineThresholdStDev(., stdev=3,thresholds = thresholds)
      # - [x] set clusters ntc control with new thresholds
      control.data.ntc %<>%
        defineClusters(., thresholds.2)
      # - [x] set file name NTC control sample 
      control.name <- paste(file.names[sample.type == "ntc" | sample.type == "neg"],"_ntc_Control",sep="",collapse="")
      output.file <- file.path(path.targets[i], paste(control.name,".png",sep=""))
      # - [x] create plot for control data
      png(filename=output.file,width = 800,height = 800)
      plot.ddpcr(x=control.data.ntc, main=control.name, max.xy=data.xy.max, thresholds=thresholds.2)
      dev.off()
      # - [x] set results <- c()
      results <- c()
      # - [x] loop through all the samples
      for(j in 1:length(files))
      {
        sample.data <- 
          read.table(file=file.path(path,files[j]),header = TRUE,sep = ",") %>%
          defineClusters(., thresholds)
        # - [x] redefine clusters with mean & stdev
        thresholds.2 <- 
          sample.data %>% 
          refineThresholdStDev(., stdev=3,thresholds = thresholds)
        # - [x] set clusters sample data with new thresholds
        sample.data %<>%
          defineClusters(., thresholds.2)
        # - [x] create plot for sample data
        output.file <- file.path(path.targets[i], paste(file.wells[j],"_",file.names[j],".png",sep=""))
        png(filename=output.file,width = 800,height = 800)
        plot.ddpcr(x=sample.data, main=paste(file.wells[j],": ",file.names[j]), max.xy = data.xy.max, thresholds = thresholds.2)
        dev.off()
        
        # - [x] Add: Well, Sample, TargetType (ch1/ch2), Target, Status concentration, 
        result <- cbind(Well=rep(as.character(file.wells[j]),2),Sample=rep(as.character(file.names[j]),2))
        result <- data.frame(result)
        result <- cbind(result, TargetType=c("Channel 1","Channel 2"))
        result <- cbind(result, Target=rep(targets[i],2))
        result <- cbind(result, Status=rep(sample.qc(x=sample.data,sample.type = sample.type[j]),2))
        result <- cbind(result, Threshold=thresholds.2)
        result <- cbind(result, calculateStatsDroplets(sample.data))
        copies.data <- CalculateStatsCopies(sample.data)
        result <- cbind(result, copies.data)
        result <- cbind(result, ngPer1ul=convertCopiesToNanogram(result$CopiesPer1ul))
        result <- cbind(result, calculateStatsRatioFract(copies.data))
        # - [ ] sample status needs to be updated
        results <- rbind(results,result)
      }
      output.file <- file.path(path.targets[i],paste(experiment,"_",targets[i], "_results.txt",sep=""))
      write.table(x = results,file = output.file,quote = FALSE,sep = "\t",row.names = FALSE)
      cat("Probe", targets[i], "has been processed.\n")
    } else {cat("Not all files for target ",targets[[i]]," are present for processing.\n", sep="")}
  }
}
####
# SCRIPT USES A DESIGN FILE TO RUN THE PIPELINE
# 
#
####

computer.name <- as.character(Sys.info()["nodename"]) 
if(computer.name == "localhost.localdomain")
{
  source.path <- "/home/dirk/Desktop/ownCloud/r_scripts/ddPCR_analysis/scripts_ddPCR"
  R.utils::sourceDirectory(source.path, modifiedOnly=FALSE);
  data.path <- "/run/user/1000/gvfs/smb-share:server=vumc-cl-fs02-01,share=microarray-faciliteit$/Data ddPCR/20160317 EGFR test sm_2016-03-17-14-59"
}
if(computer.name == "MacBook-Air-van-Dirk.local")
{ 
  source.path <- "/Users/dirkvanessen/ownCloud/r_scripts/ddPCR analysis/scripts ddPCR"
  R.utils::sourceDirectory(source.path, modifiedOnly=FALSE);
  data.path <- "/Users/dirkvanessen/ownCloud/r_scripts/ddPCR analysis/input.data/20160317 EGFR test sm_2016-03-17-14-59"
}

### create design 
analysis.path <- "/run/user/1000/gvfs/smb-share:server=vumc-cl-fs02-01,share=microarray-faciliteit$/Data ddPCR/20160317 EGFR test sm_2016-03-17-14-59"

if(file.exists(file.path(analysis.path, "design.txt")) == FALSE)
{
  createDesign(path=analysis.path)
}

test <- .rs.askForPassword(prompt = "Is design file ready (y/n)")
if(tolower(test) %in% c("y","yes","true") == TRUE)
{
  cat ("Processing ddPCR run.\n")
  ddpcr.analysis(path=analysis.path)
  
  # analysis.path <- gsub(pattern = basename(analysis.path), replacement = "", x = analysis.path)
}


  


