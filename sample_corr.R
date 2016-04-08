
source("D:\\R SCRIPTS\\ddPCR analysis\\scripts\\ddPCR.R")

path <- "W:\\Data ddPCR\\20151217 EGFR spike-in sm_2015-12-17-16-48"


    experiment  <- list.files(path, pattern = "Error.log",full.names = FALSE)
    experiment <- gsub(pattern = "Error.log", replacement = "",x = experiment)
    data.targets <- get.targets(path = path)
    # - [x] create new directory for the target
    path.targets <- create.target.folders(x = data.targets, path = path)
    targets <- names(data.targets)
    i=1
    

      sample.list  <- data.targets[[i]]
      sample.type <- get.controls(x = sample.list$Sample[duplicated(sample.list$Well)],pos = control.sample) 
      file.names <- unique(sample.list$Sample)
      file.wells <- unique(sample.list$Well)
      files <- paste(experiment,"_",file.wells,"_Amplitude.csv",sep="")
      
      positive.control.wells <- c("A01","A03","B01","B03","C01","D01","E01","F01","G01","H01")
      positive.control.wells <- unique(grep(paste(positive.control.wells,collapse="|"), files, value=TRUE))
    
        # - [x] get max Amplitude of all the files
        data.xy.max <- 
          combine.samples(path=path,files=positive.control.wells) %>%
          get.max.channels(.)
        # - [x] get positive sample and determine threshold
        control.data.pos  <- 
          positive.control.wells %>%
          combine.samples(path=path,files=.)
        # - [x] get positive control threshold
        thresholds <- 
          control.data.pos %>% 
          get.ddpcr.thresholds(., algorithm = "hist")
        # - [x] set clusters positive control with thresholds
        control.data.pos %<>%dim
          define.clusters(., thresholds)
        
        

        sample.cluster.data <- function(x, sample.size=1500, pos.percentage=50,replace = FALSE)
        {
          pos.droplets <- round(sample.size/100)*pos.percentage
          neg.droplets <- sample.size-pos.droplets
          neg.cluster <- sample(x[x[,3] == 1,1],size = neg.droplets, replace = replace)
          pos.cluster <- sample(x[x[,3] == 2,1],size = pos.droplets, replace = replace)
          cor.data <- c(neg.cluster,pos.cluster)
          cor.data <- cor.data[order(cor.data)]
          return(cor.data)
        }
       
        percentages <- c(0,10,20,30,40,50,60,70,80,90,100)
        results <- c()
        for(i in 1:length(percentages))
        {
          temp <- sample.cluster.data(x = control.data.pos, sample.size = 1500, pos.percentage = percentages[i])
          results <- cbind(results,temp)
          colnames(results)[i] <- paste("perc_",percentages[i],sep="")
        }
        
        output.file <- file.path("D:\\R SCRIPTS\\ddPCR analysis\\development","pos_sample_results.txt",sep="")
        write.table(x = results,file = output.file,quote = FALSE,sep = "\t",row.names = FALSE)
        
        
        control.data.ntc <- 
          files[sample.type == "ntc"] %>%
          combine.samples(path=path,files=.)
        
        sample <- sample(control.data.ntc[,1], 1500)
        sample <- sample[order(sample)]
        
        for(i in 1:dim(results)[2])
        {
          cat(cor(results[,i],sample),"\n")
        }
       
        
        
        # - [x] set file name control sample 
        control.name <- paste(file.names[sample.type == "pos"],"_pos_Control",sep="", collapse="")
        output.file <- file.path(path.targets[[i]], paste(control.name,".png",sep=""))
        # - [x] create plot for control data
        png(filename=output.file,width = 800,height = 800)
        plot.ddpcr(x=control.data.pos, main=control.name, max.xy=data.xy.max, thresholds=threshold)
        dev.off()

        # - [x] get ntc sample(s) 
        control.data.ntc <- 
          files[sample.type == "ntc"] %>%
          combine.samples(path=path,files=.)
        # - [ ] set clusters ntc control with thresholds
        control.data.ntc %<>%
          define.clusters(., thresholds)
        # - [x] redefine clusters with mean & stdev
        threshold.2 <- 
          control.data.ntc %>% 
          refine.clusters.stdev(., stdev=3,thresholds = thresholds)
        # - [x] set clusters ntc control with new thresholds
        control.data.ntc %<>%
          define.clusters(., thresholds.2)
        # - [x] set file name NTC control sample 
        control.name <- paste(file.names[sample.type == "ntc"],"_ntc_Control",sep="",collapse="")
        output.file <- file.path(path.targets[[i]], paste(control.name,".png",sep=""))
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
            define.clusters(., thresholds)
          # - [x] redefine clusters with mean & stdev
          thresholds.2 <- 
            sample.data %>% 
            refine.clusters.stdev(., stdev=3,thresholds = thresholds)
          # - [x] set clusters sample data with new thresholds
          sample.data %<>%
            define.clusters(., thresholds.2)
          # - [x] create plot for sample data
          sample.name <- gsub(pattern = "_Amplitude.csv",replacement="",x=files[j])
          output.file <- file.path(path.targets[[i]], paste(sample.name,".png",sep=""))
          png(filename=output.file,width = 800,height = 800)
          plot.ddpcr(x=sample.data, main=file.names[j], max.xy = data.xy.max, thresholds = thresholds.2)
          dev.off()
          
          # - [x] Add: Well, Sample, TargetType (ch1/ch2), Target, Status concentration, 
          result <- cbind(Well=rep(as.character(file.wells[j]),2),Sample=rep(as.character(file.names[j]),2))
          result <- data.frame(result)
          result <- cbind(result, TargetType=c("Channel 1","Channel 2"))
          result <- cbind(result, Target=rep(targets[i],2))
          result <- cbind(result, Status=rep(sample.qc(x=sample.data,sample.type = sample.type[j]),2))
          result <- cbind(result, Threshold=thresholds.2)
          result <- cbind(result, get.statistics.droplets(sample.data))
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
  
