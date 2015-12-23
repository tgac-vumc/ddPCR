
source("D:\\R SCRIPTS\\ddPCR analysis\\scripts\\ddPCR.R")

  ddpcr.analysis <- function(path)
  {
    experiment  <- list.files(path, pattern = "Error.log",full.names = FALSE)
    experiment <- gsub(pattern = "Error.log", replacement = "",x = experiment)
    data.targets <- get.targets(path = path)
    # - [x] create new directory for the target
    path.targets <- create.target.folders(x = data.targets, path = path)
    targets <- names(data.targets)
    
    # - [x] loop through all the targets
    for(i in 1:length(targets))
    {
      if(targets[i] == "EGFR L858R" | targets[i] == "EGFR_L858R")
      {
        control.sample <- "H1975"
      } else if(targets[i] == "EGFR T790M" | targets[i] == "EGFR_T790M")
      {
        control.sample <- "H1975"
      } else if(targets[i] == "EGFR E746_A750" | targets[i] == "EGFR_E746_A750" | targets[i] == "EGFR_E746_A750del")
      {
        control.sample <- "H1650"
      } else 
        {
          cat("No control sample has been found./n")
          break
        }

      sample.list  <- data.targets[[i]]
      sample.type <- get.controls(x = sample.list$Sample[duplicated(sample.list$Well)],pos = control.sample) 
      file.names <- unique(sample.list$Sample)
      file.wells <- unique(sample.list$Well)
      files <- paste(experiment,"_",file.wells,"_Amplitude.csv",sep="")
      
      if(sum(files %in% list.files(path)) == length(files))
        {# are all files available for analysis.
        # - [x] get max Amplitude of all the files
        data.xy.max <- 
          combine.samples(path=path,files=files) %>%
          get.max.channels(.)
        # - [x] get positive sample and determine breakpoints
        control.data.pos  <- 
          files[sample.type == "pos"] %>%
          combine.samples(path=path,files=.)
        # - [x] get positive control breakpoints
        breakpoints <- 
          control.data.pos %>% 
          get.ddpcr.breakpoints(., algorithm = "hist")
        # - [x] set clusters positive control with breakpoints
        control.data.pos %<>%
          define.clusters(., breakpoints)
        # - [x] set file name control sample 
        control.name <- paste(file.names[sample.type == "pos"],"_pos_Control",sep="", collapse="")
        output.file <- file.path(path.targets[[i]], paste(control.name,".png",sep=""))
        # - [x] create plot for control data
        png(filename=output.file,width = 800,height = 800)
        plot.ddpcr(x=control.data.pos, main=control.name, max.xy=data.xy.max, breakpoints=breakpoints)
        dev.off()

        # - [x] get ntc sample(s) 
        control.data.ntc <- 
          files[sample.type == "ntc"] %>%
          combine.samples(path=path,files=.)
        # - [ ] set clusters ntc control with breakpoints
        control.data.ntc %<>%
          define.clusters(., breakpoints)
        # - [x] redefine clusters with mean & stdev
        breakpoints.2 <- 
          control.data.ntc %>% 
          refine.clusters.stdev(., stdev=3,breakpoints = breakpoints)
        # - [x] set clusters ntc control with new breakpoints
        control.data.ntc %<>%
          define.clusters(., breakpoints.2)
        # - [x] set file name NTC control sample 
        control.name <- paste(file.names[sample.type == "ntc"],"_ntc_Control",sep="",collapse="")
        output.file <- file.path(path.targets[[i]], paste(control.name,".png",sep=""))
        # - [x] create plot for control data
        png(filename=output.file,width = 800,height = 800)
        plot.ddpcr(x=control.data.ntc, main=control.name, max.xy=data.xy.max, breakpoints=breakpoints)
        dev.off()
        # - [x] set results <- c()
        results <- c()
        # - [x] loop through all the samples
        for(j in 1:length(files))
        {
          sample.data <- 
            read.table(file=file.path(path,files[j]),header = TRUE,sep = ",") %>%
            define.clusters(., breakpoints)
          # - [x] redefine clusters with mean & stdev
          breakpoints.2 <- 
            sample.data %>% 
            refine.clusters.stdev(., stdev=3,breakpoints = breakpoints)
          # - [x] set clusters sample data with new breakpoints
          sample.data %<>%
            define.clusters(., breakpoints.2)
          # - [x] create plot for sample data
          sample.name <- gsub(pattern = "_Amplitude.csv",replacement="",x=files[j])
          output.file <- file.path(path.targets[[i]], paste(sample.name,".png",sep=""))
          png(filename=output.file,width = 800,height = 800)
          plot.ddpcr(x=sample.data, main=file.names[j], max.xy = data.xy.max, breakpoints = breakpoints.2)
          dev.off()
          
          # - [x] Add: Well, Sample, TargetType (ch1/ch2), Target, Status concentration, 
          result <- cbind(Well=rep(as.character(file.wells[j]),2),Sample=rep(as.character(file.names[j]),2))
          result <- data.frame(result)
          result <- cbind(result, TargetType=c("Channel 1","Channel 2"))
          result <- cbind(result, Target=rep(targets[i],2))
          result <- cbind(result, Status=rep(sample.qc(x=sample.data,sample.type = sample.type[j]),2))
          result <- cbind(result, Threshold=breakpoints.2)
          result <- cbind(result, get.statistics.droplets(sample.data))
            # - [x] colnames of droplet count data is changed after data.frame conversion
          copies.data <- get.statistics.copies(sample.data)
          result <- cbind(result, copies.data)
          result <- cbind(result, ngPer1ul=convert.copies.to.ng(result$CopiesPer1ul))
          result <- cbind(result, get.statistics.ratio.fract(copies.data))
          # - [ ] sample status needs to be updated
          results <- rbind(results,result)
        }
        output.file <- file.path(path.targets[[i]],paste(experiment,"_",targets[[i]], "_results.txt",sep=""))
        write.table(x = results,file = output.file,quote = FALSE,sep = "\t",row.names = FALSE)
        cat("Probe", targets[i], "has been processed.\n")
        } else {cat("Not all files for target ",targets[[i]]," are present for processing.\n", sep="")}
    }
  }
  analysis.path <- choose.dir()
  ddpcr.analysis(path=analysis.path)
  
