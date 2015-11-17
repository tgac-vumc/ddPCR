input.path <- "D:\\R SCRIPTS\\ddPCR analysis" #work
#input.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis"
folders <- c("archive","input.data","output.data","output.plot","scripts","scripts.log")

library(magrittr)
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

archive <- function()
{
  # "design file" with sample information
  # name  file/well type  probe FDR_positive  FDR_rain
  # name = sample name
  # name/well = sample well name /sample file nae
  # sample type = 'positive', 'sample', 'negative'
  # probe = probe name to select the 'probe file
  # - [x] get cluster data - mean and sd for the clusters 

  define.color <- function(x, breakpoint.ch1, breakpoint.ch2,density=40)
  {
    # - [x] find cluster notation BioRad 
    ddpcr.colors <- paste(c("#000000","#0033FF","#FF6600","#00CC00"), as.character(density), sep="")
    results <- rep(NA,dim(x)[1])
    results[x[,1] < breakpoint.ch1 & x[,2] < breakpoint.ch2] <- ddpcr.colors[1] # ch1-ch2- : cluster 1
    results[x[,1] > breakpoint.ch1 & x[,2] < breakpoint.ch2] <- ddpcr.colors[2] # ch1+ch2- : cluster 2
    results[x[,1] > breakpoint.ch1 & x[,2] > breakpoint.ch2] <- ddpcr.colors[3] # ch1+ch2+ : cluster 3
    results[x[,1] < breakpoint.ch1 & x[,2] > breakpoint.ch2] <- ddpcr.colors[4] # ch1-ch2+ : cluster 4
    return(results)
  }
  clusters.mean.sd <- function(x,na.rm=TRUE,breakpoint)
  {
    clusters <- c(mean(x[x < breakpoint],na.rm=na.rm), sd(x[x < breakpoint],na.rm=na.rm))
    clusters <- rbind(clusters,c(mean(x[x > breakpoint],na.rm=na.rm), sd(x[x > breakpoint],na.rm=na.rm)))
    clusters <- cbind(clusters,clusters[,2]*3)
    rownames(clusters) <- c("cluster1","cluster2");colnames(clusters) <- c("mean","sd","3*sd")
    return(clusters)
  }
  mean.sd.ch1 <- clusters.mean.sd(x = control.data[,1],breakpoint = breakpoint.ch1)
  mean.sd.ch2 <- clusters.mean.sd(x = control.data[,2],breakpoint = breakpoint.ch2)
  
  # - [x] find the cut-off values for negative cluster
  times.sd <- 4
  
  ch1.neg.cutoff <- mean.sd.ch1[1,1]+(mean.sd.ch1[1,2] * times.sd)
  abline(h=ch1.neg.cutoff,col="black")
  ch2.neg.cutoff <- mean.sd.ch2[1,1]+(mean.sd.ch2[1,2]* times.sd)
  abline(v=ch2.neg.cutoff,col="black")
  
  # - [x] find the cut-off values for ch1 positive
  ch1.neg.cutoff <- mean.sd.ch1[2,1]-(mean.sd.ch1[2,2] * times.sd)
  abline(h=ch1.neg.cutoff,col="blue")
  ch2.neg.cutoff <- mean.sd.ch2[1,1]+(mean.sd.ch2[1,2]* times.sd)
  abline(v=ch2.neg.cutoff,col="blue")
}

ddpcr.pipeline <- function()
  {
  # - [ ] PIPELINE FUNCTIONS
  mgsub <- function(pattern, replacement, x, ...) 
  {
    if (length(pattern)!=length(replacement)) {
      stop("pattern and replacement do not have the same length.")
    }
    result <- x
    for (i in 1:length(pattern)) {
      result <- gsub(pattern[i], replacement[i], result, ...)
    }
    result
  }
  create.design.file <- function(path,probe="probe")
  {
    # - [x] create design file for experiment
    amplitude.files <- list.files(path = path, pattern = "_Amplitude.csv")
    sample.names <- gsub(pattern = "_Amplitude.csv",x = amplitude.files, replacement = "")
    design <- matrix(data = "", nrow = length(amplitude.files), ncol = 4)
    colnames(design) <- c("Name","File","Type","Probe")
    design[,1] <- sample.names
    design[,2] <- amplitude.files
    design[,3] <- "pos-neg-sample"
    design[,4] <- probe
    output.file <- file.path(path,"design.txt")
    write.table(file = output.file,x = design,quote = FALSE,sep = "\t",row.names = FALSE)
  }
  add.probe.data <- function(path,Name,Date,BreakPoint1,BreakPoint2,PosClusterAmp,NegClusterAmp,PosClusterSD,NegClusterSD,PosFDR,RainFDR)
  {
    Name <- tolower(Name)
    probe.file <- file.path(path,paste(name,".Rdata",sep=""))
    if(!file.exists(probe.file)){
      col.names <- c("Ch1","Ch2")
      probe.data <- list( Name = Name,
                          Date = Date,
                          BreakPoint1 = matrix(data = BreakPoint1, nrow = 1, ncol = 2,dimnames = list(NULL,col.names)),
                          BreakPoint2 = matrix(data = BreakPoint2, nrow = 1, ncol = 2,dimnames = list(NULL,col.names)),
                          PosClusterAmp = matrix(data = PosClusterAmp, nrow = 1, ncol = 2,dimnames = list(NULL,col.names)),
                          NegClusterAmp = matrix(data = NegClusterAmp, nrow = 1, ncol = 2,dimnames = list(NULL,col.names)),
                          PosClusterSD = matrix(data = PosClusterSD, nrow = 1, ncol = 2,dimnames = list(NULL,col.names)),
                          NegClusterSD = matrix(data = NegClusterSD, nrow = 1, ncol = 2,dimnames = list(NULL,col.names)),
                          PosFDR = matrix(data = PosFDR, nrow = 1, ncol = 2,dimnames = list(NULL,col.names)),
                          RainFDR = matrix(data = RainFDR, nrow = 1, ncol = 2,dimnames = list(NULL,col.names))
                          )
      save(probe.data, file=probe.file)
      # - [x] create basic probe file
      # - [ ] save probe file
    }else{
      load(probe.file)
        probe.data$Name <- rbind(probe.data$Name,Name)
        probe.data$Date <- rbind(probe.data$Date,Date)
        probe.data$BreakPoint1 <- rbind(probe.data$BreakPoint1,BreakPoint1)
        probe.data$BreakPoint2 <- rbind(probe.data$BreakPoint2,BreakPoint2)
        probe.data$PosClusterAmp <- rbind(probe.data$PosClusterAmp,PosClusterAmp)
        probe.data$NegClusterAmp <- rbind(probe.data$NegClusterAmp,NegClusterAmp)
        probe.data$PosClusterSD <- rbind(probe.data$PosClusterSD,PosClusterSD)
        probe.data$NegClusterSD <- rbind(probe.data$NegClusterSD,NegClusterSD)
        probe.data$PosFDR <- rbind(probe.data$PosFDR,PosFDR)
        probe.data$RainFDR <- rbind(probe.data$RainFDR,RainFDR)
        save(probe.data, file=probe.file)
      # - [ ] read probe file
      # - [ ] add probe data to file
      # - [ ] save probe file
    }
  }
  add.probe.data.v2 <- function(path,name=NA,date=NA,breakpoints=c(NA,NA))
  {
    name <- tolower(name)
    probe.file <- file.path(path,paste(name,".Rdata",sep=""))
    if(!file.exists(probe.file)){
      col.names <- c("Ch1","Ch2")
      probe.data <- list( name = name,
                          date = date,
                          breakpoints = matrix(data = breakpoints, nrow = 1, ncol = 2,dimnames = list(NULL,col.names))
      )
      save(probe.data, file=probe.file)
    }else{
      load(probe.file)
      probe.data$name <- rbind(probe.data$name,name)
      probe.data$date <- rbind(probe.data$date,date)
      probe.data$breakpoints <- rbind(probe.data$breakpoints,breakpoints)
      save(probe.data, file=probe.file)
    }
    return(probe.data)
  }
  read.design.file <- function(path,pattern=NULL)
  {
    design.file <- list.files(path=path,pattern=pattern)
    if(length(design.file) > 1){
      stop("multiple design files found. Unable to continue.\n"
      )
    }else{
      design <- read.table(file = file.path(path,design.file),header = TRUE,sep = "\t")
      return(design)
    }
  }
  combine.samples <- function(path,files)
  {
    combined.data <- NULL
    for(i in 1:length(files))
    { 
      sample.data <- read.table(file=file.path(path,files[i]),header = TRUE,sep = ",")
      combined.data <- rbind(combined.data,sample.data )
    }
    return(combined.data)
  }
  get.breakpoint.kmeans <- function(x,nClusters=2)
  { # use kmeans function
    x <- as.numeric(x)
    result <- NULL
    breakpoint <- kmeans(x=x,centers=nClusters)$centers
    if(dim(breakpoint)[1] == 2){result <- mean(breakpoint)}
    return(result)
  }
  get.breakpoint.ranges <- function(x)
  { # use min/max function
    x <- as.numeric(x)
    result <- NULL
    result <-  (max(x) - min(x))/2
    return(result)
  }
  get.breakpoint.hist <- function(x)
  { # use hist function to determine locations with no droplets (middle?)
    x <- as.numeric(x)
    hist.data <- rbind(hist(x, breaks=15, plot=FALSE)$mids, hist(x, breaks=15, plot=FALSE)$counts)
    hist.data <- hist.data[,-c(1:2,14:16)]
    result <- mean(hist.data[1,hist.data[2,] == min(hist.data[2,])])
    return(result)
  }
  get.ddpcr.breakpoints.kmeans <- function(x)
  { # use kmeans on x & y together, and set for finding 4 clusters
    results <- kmeans(x[,1:2], 4)
    return(results)
  }
  get.ddpcr.breakpoints <- function(x, algorithm = "hist")
  { 
    if(tolower(algorithm) == "hist" | tolower(algorithm) == "histogram")
    {
      results <- c(breakpoint.ch1 = get.breakpoint.hist(x = x[,1]) ,breakpoint.ch2 = get.breakpoint.hist(x = x[,2]))
    }
    if(tolower(algorithm) == "ranges")
    {
      results <- c(breakpoint.ch1 = get.breakpoint.ranges(x = x[,1]) ,breakpoint.ch2 = get.breakpoint.ranges(x = x[,2]))
    }
    if(tolower(algorithm) == "kmeans")
    {
      results <- c(breakpoint.ch1 = get.breakpoint.kmeans(x = x[,1]) ,breakpoint.ch2 = get.breakpoint.kmeans(x = x[,2]))
    }
    
    return(results)
  }
  define.clusters <- function(x, breakpoints)
  {
    if (length(breakpoints) != 2) {
      stop("breakpoints must have a length of 2.\n")
    }
    # - [x] find cluster notation BioRad 
    results <- rep(NA,dim(x)[1])
    results[x[,1] < breakpoints[1] & x[,2] < breakpoints[2]] <- 1 # ch1-ch2- : cluster 1
    results[x[,1] > breakpoints[1] & x[,2] < breakpoints[2]] <- 2 # ch1+ch2- : cluster 2
    results[x[,1] > breakpoints[1] & x[,2] > breakpoints[2]] <- 3 # ch1+ch2+ : cluster 3
    results[x[,1] < breakpoints[1] & x[,2] > breakpoints[2]] <- 4 # ch1-ch2+ : cluster 4
    x[,3] <- results
    return(x)
  }
  define.color <- function(x,density=40)
  {
    # - [x] find cluster notation BioRad 
    ddpcr.colors <- paste(c("#000000","#FF6600","#00CC00","#0033FF"), as.character(density), sep="")
    x <- as.character(x)
    x <- mgsub(pattern = c("1","3","4","2"),replacement = ddpcr.colors,x=x)
    return(x)
  }
  dropletcount.clusters <- function(x)
  {
    results <- NULL
    results <- list(clusters=c(cluster.1=sum(x == 1),cluster.2=sum(x == 2),cluster.3=sum(x == 3),cluster.4=sum(x == 4)))
    results  <-paste("Ch1-Ch2-:",results$clusters[1],
                             "   Ch1+Ch2-:",results$clusters[2],
                             "   Ch1+Ch2+:",results$clusters[3],
                             "   Ch1-Ch2+:",results$clusters[4], sep="")
    return(results)
  }
  get.max.channels <- function(x)
  {
    results <- c(Ch1.max = round(max(x[,1])+100) ,Ch2.max = round(max(x[,2])+100))
    return(results)
  }
  plot.ddpcr <- function(x,dotres=0.7,main="ddPCR",pch=16,colors="ddpcr",density=60,breakpoints=NULL,max.xy=NULL,verbose=FALSE)
  {
    if(length(max.xy) != 2) {
      xmax <- max(x[,2])
      ymax <- max(x[,1])
    } else {
      xmax <- max.xy[2]
      ymax <- max.xy[1]
    }
    col.vec <- define.color(x = x[,3],density = density)
    plot(y=x[,1],x=x[,2], cex=dotres, col=col.vec, ylab="Ch1 Amplitude",xlab="Ch2 Amplitude", pch=pch, main=main,
         xlim=c(0,xmax),ylim=c(0,ymax))
    sub.text <- dropletcount.clusters(x=x[,3])
    mtext(side = 3,text = sub.text, cex = 0.8)
    if (length(breakpoints) != 2) {
      if(verbose == TRUE){cat("No breakpoint data has been given. Data will not be plotted.")}
    }else{
      abline(h=breakpoints[1], col="red") # channel 1
      abline(v=breakpoints[2], col="red") # channel 2
    }
  }
  
  
  
  ddpcr.analysis <- function(path,probe.path)
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
      dev.off()
    }
  }
  # - [ ] PIPELINE SETUP
  create.design.file(project.path,probe = "test")
  
  ##### SET PROJECT PATH
  probe.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/probe.data"
  
  project.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\E746_A750del"
  ddpcr.analysis(path = project.path)
  project.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\L858R"
  ddpcr.analysis(path = project.path)
  project.path <- "D:\\R SCRIPTS\\ddPCR analysis\\input.data\\T790M"
  ddpcr.analysis(path = project.path)
  
  project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/L858R"
  ddpcr.analysis(path = project.path, probe.path = probe.path)
  project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/E746_A750del"
  ddpcr.analysis(path = project.path, probe.path = probe.path)
  project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/T790M"
  ddpcr.analysis(path = project.path, probe.path = probe.path)
  project.path <- "/Users/dirkvanessen/Desktop/ddPCR analysis/input.data/test"
  ddpcr.analysis(path = project.path, probe.path = probe.path)
  
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
 }
  
  

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
  
  

