input.path <- "D:\\R SCRIPTS\\ddPCR analysis" # work
#input.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis" # home
folders <- c("archive","input.data","output.data","output.plot","scripts","scripts.log")
# FUNCTIONS
get_comments = function(filename){
  is_assign = function(expr) as.character(expr) %in% c("<-", "<<-", "=", "assign")
  is_function = function(expr) is.call(expr) && is_assign(expr[[1L]]) && is.call(expr[[3L]]) && expr[[3L]][[1L]] == quote(`function`)
  src = parse(filename, keep.source = TRUE)
  functions = Filter(is_function, src)
  fun_names = as.character(lapply(functions, `[[`, 2L))
  # - [x] extract all comments
  r = setNames(lapply(attr(functions, "srcref"), grep, pattern = "^\\s*#", value = TRUE), fun_names)
  # - [x] remove leading spaces and comment sign '#'
  r = lapply(r, function(x) sub(pattern = "^\\s*#", replacement = "", x = x))
  # - [x] keep only markdown checkboxes like " - [ ] " or " - [x] "
  r = lapply(r, function(x) x[nchar(x) >= 7L & substr(x, 1L, 7L) %in% c(" - [ ] "," - [x] ")])
  # - [x] return only non empty results
  r[as.logical(sapply(r, length))]
}
make_doc = function(path = "R", files, package, dest){
  if(!missing(package)) path = system.file(path, package=package)
  stopifnot(file.exists(path))
  if(missing(files)) files = list.files(path, pattern = "\\.R$")
  if(!length(files)){
    warning(paste0("No files to process in ",path,"."))
    return(invisible())
  }
  if(!all(sapply(file.path(path, files), file.exists))) stop(paste0("Processing stopped as some files not exists: ", paste(files[!sapply(file.path(path, files), file.exists)], collapse=", "),"."))
  r = setNames(lapply(file.path(path, files), get_comments), files)
  r = r[as.logical(sapply(r, length))]
  if(missing(dest)) return(r)
  if(!file.exists(dirname(dest))) dir.create(dirname(dest), recursive=TRUE)
  if(file.exists(dest)) file.rename(dest, paste0(dest,"_backup"))
  invisible(lapply(names(r), function(filename){
    cat(c("",paste("###", filename)), sep = "\n", file = dest, append = file.exists(dest))
    lapply(names(r[[filename]]), function(funname){
      cat(c("",paste("####", funname),""), sep = "\n", file = dest, append = TRUE)
      cat(r[[filename]][[funname]], sep = "\n", file = dest, append = TRUE)
    })
  }))
  if(file.exists(paste0(dest,"_backup"))) file.remove(paste0(dest,"_backup"))
  invisible(dest)
}
create.folders <- function(path,folders){
  if(!file.exists(path))
  {
    stop ("Input folder does not exist. \n")
  }else
  {
    for(i in 1:length(folders))
    {
      if(!file.exists(file.path(path,as.character(folders[i]))))
      {
        dir.create(file.path(path,as.character(folders[i])))
      }
    }
  }
}
set.paths <- function(path="", folders){
  if(path == "" | class(path) == "numeric"){ stop
  } else {
    folders <- c("original",folders)
    paths <- NULL
    for(i in 1:length(folders))
    {
      if(i == 1)
      {
        paths <- list(file.path(path))
      }
      if(i > 1)
      {
        paths <- c(paths,file.path(path,as.character(folders[i])))
      }
    }
    names(paths) <- as.character(folders)
    return(paths)
  }
}
# RUN FUNCTIONS
create.folders(path=input.path,folders=folders)
path <- set.paths(path=input.path, folders=folders)
log.file <- paste(format(Sys.time(), "%Y%m%d-%H%M"),"_script_checklist.md",sep="")
make_doc(path=path$scripts,dest = file.path(path$scripts.log,log.file))
# END

  source(file.path(path$scripts,"ddPCR.R"))
# - [ ] create project design
  create.design.file(project.path,probe = "test")
  
# - [ ] PIPELINE SETUP
  ddpcr.analysis <- function(path,project.path)
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
    col.vec <- 
      control.data$Cluster %>% 
      define.color(., density=60)
    # - [x] add probe data
    # all.probe.data <- add.probe.data.v2(path = probe.path, name = exp.design[1,4], date = format(Sys.time(), "%Y-%m-%d"), breakpoints = breakpoints)
   
    # - [x] get text for plotting
    droplet.count <- 
      control.data$Cluster %>%
      dropletcount.clusters(.)
    # - [x] set file name control sample 
    control.name <- paste(strsplit(x = as.character(exp.design$File[1]),split = "_")[[1]][1],"_Controls",sep="")
    output.file <- file.path(project.path, paste(control.name,".png",sep=""))
    # - [x] create plot for control data
    png(filename=output.file,width = 800,height = 800)
    plot.ddpcr(x=control.data, main=control.name, max.xy=data.xy.max, breakpoints=breakpoints)
    dev.off()
    
    # - [x] retrieve sample files (all files)
    sample.files <- exp.design[,2]
    # - [ ] analyse sample files
    for(i in 1:length(sample.files))
    {
      sample.data <- 
        read.table(file=file.path(path,sample.files[i]),header = TRUE,sep = ",") %>%
        define.clusters(., breakpoints)
      col.vec <- 
        sample.data$Cluster %>% 
        define.color(., density=60)
      droplet.count <- 
        sample.data$Cluster %>%
        dropletcount.clusters(.)
      sample.name <- gsub(pattern = "_Amplitude.csv",replacement="",x=sample.files[i])
      output.file <- file.path(project.path, paste(sample.name,".png",sep=""))
      png(filename=output.file,width = 800,height = 800)
      plot.ddpcr(x=sample.data, main=sample.name, max.xy = data.xy.max, breakpoints = breakpoints)
      #plot.cutoffs(mean.4sd.cutoff(sample.data))
      dev.off()
    }
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
    ddpcr.analysis(path = project.path, probe.path = probe.path)
    project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/E746_A750del"
    ddpcr.analysis(path = project.path, probe.path = probe.path)
    project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/T790M"
    ddpcr.analysis(path = project.path, probe.path = probe.path)
    project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/test"
    ddpcr.analysis(path = project.path, probe.path = probe.path)
  }
  
  
  new.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\20151125_790"
  ddpcr.analysis(path = new.path, project.path=new.path)
  new.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\20151125_746_del"
  ddpcr.analysis(path = new.path, project.path=new.path)
  
 get.statistics <- function(){
   new.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\20151125_858"
   exp.design <- read.design.file(path=new.path,pattern="design")
   sample.files <- exp.design[,2]
   # - [ ] analyse sample files
   
  i = 1 # sample selection
  breakpoints <- c(2139,953)
     sample.data <- 
       read.table(file=file.path(new.path,sample.files[i]),header = TRUE,sep = ",") %>%
       define.clusters(., breakpoints)
  
  cluster.1 <- sum(sample.data[,3] %in% 1)
  cluster.2 <- sum(sample.data[,3] %in% 2)
  cluster.3 <- sum(sample.data[,3] %in% 3)
  cluster.4 <- sum(sample.data[,3] %in% 4)
  pos.ch1 <- sum(cluster.2, cluster.3)
  pos.ch2 <- sum(cluster.4, cluster.3)
  total.count <- sum(cluster.1, cluster.2, cluster.3, cluster.4)
  
  result <- concentration(negCount = (total.count-pos.ch1),Count = total.count) # channel 1
 c(result, (result + (moments(c(pos.ch1,total.count))[3]*1)))
  
  result <- concentration(negCount = (total.count-pos.ch2),Count = total.count) # channel 2
  c(result, (result + (moments(c(pos.ch2,total.count))[3]*1)))
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
  
  

