input.path <- "D:\\R SCRIPTS\\ddPCR analysis" # work
#input.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis" # home
folders <- c("archive","input.data","output.data","output.plot","scripts","scripts.log")

# FUNCTIONS
computer.name <- as.character(Sys.info()["nodename"]) 
if(computer.name == "PA-PC807")
{
  source("C:\\Documents and Settings\\h.vanessen\\ownCloud\\R_my_functions\\_R_my_functions.R")
}
if(computer.name == "iMAC")
{ 
  source("/ownCloud/R_my_functions/_R_my_functions.R")
}

# RUN FUNCTIONS
create.folders(path=input.path,folders=folders)
path <- set.paths(path=input.path, folders=folders)
log.file <- paste(format(Sys.time(), "%Y%m%d-%H%M"),"_script_checklist.md",sep="")
make_doc(path=path$scripts,dest = file.path(path$scripts.log,log.file))
# END

  source(file.path(path$scripts,"ddPCR.R"))

# - [ ] PIPELINE SETUP
  ddpcr.analysis <- function(path)
  {
    # - [x] start global analysis for experiment
    exp.design <- read.design.file(path=path,pattern="design")
    data.xy.max <- 
      combine.samples(path=path,files=exp.design$File) %>%
      get.max.channels(.)
    # - [x] analyze control files
    control.data <- 
      exp.design[exp.design$Type == "pos" | exp.design$Type == "neg",2] %>%
      combine.samples(path=path,files=.)
    breakpoints <- 
      control.data %>% 
      get.ddpcr.breakpoints(., algorithm = "hist")
    control.data %<>%
      define.clusters(., breakpoints)
    # - [x] set file name control sample 
    control.name <- paste(strsplit(x = as.character(exp.design$File[1]),split = "_")[[1]][1],"_Controls",sep="")
    output.file <- file.path(path, paste(control.name,".png",sep=""))
    # - [x] create plot for control data
    png(filename=output.file,width = 800,height = 800)
    plot.ddpcr(x=control.data, main=control.name, max.xy=data.xy.max, breakpoints=breakpoints)
    dev.off()
    
    # - [x] retrieve sample files (all files)
    sample.files <- exp.design[,2]
    results <- NULL
    
    # - [ ] analyse sample files
    for(i in 1:length(sample.files))
    {
      sample.data <- 
        read.table(file=file.path(path,sample.files[i]),header = TRUE,sep = ",") %>%
        define.clusters(., breakpoints)
      sample.name <- gsub(pattern = "_Amplitude.csv",replacement="",x=sample.files[i])
      output.file <- file.path(path, paste(sample.name,".png",sep=""))
      png(filename=output.file,width = 800,height = 800)
      plot.ddpcr(x=sample.data, main=sample.name, max.xy = data.xy.max, breakpoints = breakpoints)
      dev.off()
      # - [ ] sample names is not correct
      results <- rbind(results, get.statistics(x = sample.data,sample = exp.design[i,1],input.file = exp.design[i,2],target = exp.design[i,4],breakpoints = breakpoints))
    }
    output.file <- file.path(path,paste("ddPCR_results_" ,basename(path),".txt",sep=""))
    write.table(x = results,file = output.file,quote = FALSE,sep = "\t",row.names = FALSE)
    return(results)
  }
  
  ##### SET PROJECT PATH
  path.work <- function()
  {
    probe.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/probe.data"
    project.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\E746_A750del"
    ddpcr.analysis(path = project.path)
    project.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\L858R"
    ddpcr.analysis(path = project.path)
    project.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\T790M"
    ddpcr.analysis(path = project.path)
 } 
  path.home <- function()
  {
    project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/L858R"
    ddpcr.analysis(path = project.path)
    project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/E746_A750del"
    ddpcr.analysis(path = project.path, probe.path = probe.path)
    project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/T790M"
    ddpcr.analysis(path = project.path, probe.path = probe.path)
    project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/test"
    ddpcr.analysis(path = project.path, probe.path = probe.path)
  }
  
  new.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\spike-in L858R"
  ddpcr.analysis(path = new.path)
  new.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\spike-in T790M"
  ddpcr.analysis(path = new.path)

  
  ###################
  input.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\20151202 EGFR spike-in sm_2015-12-02-15-25"
  path <- input.path
  ddpcr.analysis.v2 <- function(path)
    {
    
    experiment  <- list.files(path, pattern = "Error.log",full.names = FALSE)
    experiment <- gsub(pattern = "Error.log", replacement = "",x = experiment)
    data.targets <- get.targets(path = path)
    # - [x] create new directory for the target
    path.targets <- create.target.folders(x = data.targets, path = path)
    targets <- names(data.targets)
    control.sample <- "H1975"
    
    for(i in 1:length(targets))
    {
      sample.list  <- data.targets[[i]]
      sample.type <- get.controls(x = sample.list$Sample[duplicated(sample.list$Well)],pos = "H1975") 
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
          control.data %>% 
          get.ddpcr.breakpoints(., algorithm = "hist")
        # - [x] set clusters positive control with breakpoints
        control.data.pos %<>%
          define.clusters(., breakpoints)
        # - [x] set file name control sample 
        control.name <- paste(file.names[sample.type == "pos"],"_pos_Control",sep="")
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
        # - [x] set file name NTC control sample 
        control.name <- paste(file.names[sample.type == "ntc"],"_ntc_Control",sep="")
        output.file <- file.path(path.targets[[i]], paste(control.name,".png",sep=""))
        # - [x] create plot for control data
        png(filename=output.file,width = 800,height = 800)
        plot.ddpcr(x=control.data.ntc, main=control.name, max.xy=data.xy.max, breakpoints=breakpoints)
        dev.off()
        # - [ ] set results <- c()
        
        # - [ ] create loop for analysis of all the samples 
          # - [ ] if ntc sample, check for positive droplets, , else check
          # - [ ] if pos sample, check for droplets in all clusters, else check
          # - [ ] create plot for sample 
        
          # - [ ] get statistics of each sample (cbind)
            # - [ ] 
      
        
        # - [ ]
       
        }
    
    
    data.xy.max <- 
      combine.samples(path=input.path,files=filenames) %>%
      get.max.channels(.)
    
    # get file names
    data.xy.max <-
    combine.samples(path=path,files=exp.design$File) %>%
      get.max.channels(.)
    
    }

}
 

 
 

  
  ### start CONTROL ANALYSIS for analysis
  # - [ ] change function 'add.probe.data'
  # - [ ] create check if pos and neg control are available
  # - [ ] retrieve control data
  # - [ ] compare control data
  # - [ ] add control data
  # - [ ] save control data
  
  ### STUFF STILL TO DO
  # - [ ] save plots in plot folder automatically
  # - [ ] save processed data files in folder
  # - [ ] change breakpoint v4: use mean neg droplets as minimum value for finding breakpoint
  # - [ ] use hist & hist$density to determine the cluster mean of the negative cluster. With the breakpoints, find the negative cluster 3sd values.
 
  
  

  # probe file with probe information
  # date of each data set
  # breakpoint channel 1, multiple (learning)
  # breakpoint channel 2, multiple (learning)
  # cluster amplitude channel 1, multiple (learning)
  # cluster amplitude channel 2, multiple (learning)
  # cluster sd channel 1, multiple (learning)
  # cluster sd channel 2, multiple (learning)
  # FDR_positive = the number of positive droplets in the negative control in the 'positive cluster'
  # FDR_rain = the number of droplets in the negative control in the 'rain'
  
  # - [ ] ideas
  # read all data from all wells and find the breakpoint for all (both channels)
  # there must be 4 clusters for all the data (positive/negative/samples)
  # if auto clustering is done (without set breakpoint), is the breakpoint then at the same level?
  # is the difference between positive amplitude and negative amplitude the same for the control?
  
  # - [ ] function: create basic design file based on files in a folder
  # strsplit and find if they are similar...
  # combine data from multiple plates option?
  # - [ ] function: read design
  # - [ ] function: check files needed for analysis : sample files, probe file
  # - [ ] function: read sample files
  # - [ ] function: check positive / negative controls available
  # - [ ] function: compare positive / negative controls to probe file data
  # - [ ] function: normalize the data based on negative probes ? mean/median/mode option?
  # - [ ] function: calculate calls
  # - [ ] function: parse probe controls and store data in probe data file
  # - [ ] find the min & max for all droplets 
  # use min / max values for plotting of all the wells in this analysis set
  
  # - [ ] combine positive and negative control 
  # calculate calls (see below)
  # calculate the breakpoint
  # calculate the amplitude of cluster 1 & 2, channel 1
  # calculate the amplitude of cluster 1 & 2, channel 2
  # are there 3 or 4 clusters for the positive/negative controls
  # should all samples be normalized based on the mean/median of the negative cluster data?
  # should a flag class be used based on a 2-D plot? 
  # Ch1+Ch2+	Ch1+Ch2-	Ch1-Ch2+	Ch1-Ch2-  Ch1+Ch2R Ch1+CH2R
  
  # negative droplet channel 1 & negative channel 2
  # negative droplet channel 1 & rain channel 2
  # negative droplet channel 2 & rain channel 1
  # negative droplet channel 1 & positive droplet channel 2
  # negative droplet channel 2 & positive droplet channel 1
  # positive droplet channel 1 & rain channel 2
  # positive droplet channel 2 & rain channel 1
  # positive droplet channel 1 & positive droplet channel 2 
  # rain channel 1 & rain channel 2
  
  # - [ ] for each sample
  # find the clusters in the sample based on the positive and negative cluster data and the breakpoint
  # is the positive cluster amplitude at a similar height as the control?
  # calculate the amount of positive droplets -> 1
  # calculate the amount of negative droplets -> -1
  # calculate the amount of droplets rain -> 0
  # add a column to flag the data with the calls (-1,0,1)
  
  # - [ ] calculate the FDR base on the positive and negative controls
  # add calculations/results to probe file 
  # how many of these results should be stored
  
  # - [ ] FDR ANALYSIS
  # read in all positive and negative controls
  # calculate the FDR for all negative samples
  
  # - [ ] 
  # - [ ] 
  # - [ ] 
  
  

