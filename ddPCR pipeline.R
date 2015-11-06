input.path <- "D:\\R SCRIPTS\\ddPCR analysis" #work
project <- "Test_EGFR_2015-11-03-16-39"
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
  # - [x] set path
  # - [x] set folders
  # - [x] for loooooop
  for(i in 1:length(folders))
  {
    if(!file.exists(file.path(path,as.character(folders[i]))))
    {
      dir.create(file.path(path,as.character(folders[i])))
    }
  }
  
  # - [x] check file.exists
  # - [x] create folder
  # - [x] does it work?
}
set.paths <- function(path=""){
  if(path == "" | class(path) == "numeric"){ stop
  } else {
    # - [x] set paths
    paths <- list(archive=file.path(path,"archive"),
                  input.data=file.path(path,"input.data"),
                  output.data=file.path(path,"output.data"),
                  output.plot=file.path(path,"output.plot"),
                  scripts=file.path(path,"scripts"),
                  scripts.log= file.path(path,"scripts.log"))
    # - [x] make list
    return(paths)
    # - [x] return list
  }
}
# RUN FUNCTIONS
create.folders(path=input.path,c("archive","input.data","output.data","output.plot","scripts","scripts.log"))
path <- set.paths(input.path)
log.file <- paste(format(Sys.time(), "%Y%m%d-%H%M"),"_script_checklist.md",sep="")
make_doc(path=path$scripts,dest = file.path(path$scripts.log,log.file))


ddpcr.pipeline <- function()
  {
  
  # - [ ] create files that are needed for the pipeline

  create.design.file <- function(path)
    {
    # - [x] create design file for experiment
    amplitude.files <- list.files(path = path, pattern = "_Amplitude.csv")
    sample.names <- gsub(pattern = "_Amplitude.csv",x = amplitude.files, replacement = "")
    design <- matrix(data = "", nrow = length(amplitude.files), ncol = 4)
    colnames(design) <- c("Name","File","Type","Probe")
    design[,1] <- sample.names
    design[,2] <- amplitude.files
    design[,3] <- "pos-neg-sample"
    design[,4] <- "probe_name"
    output.file <- file.path(path,"design.txt")
    write.table(file = output.file,x = design,quote = FALSE,sep = "\t",row.names = FALSE)
  }
  # "design file" with sample information
  # name  file/well type  probe FDR_positive  FDR_rain
  # name = sample name
  # name/well = sample well name /sample file nae
  # sample type = 'positive', 'sample', 'negative'
  # probe = probe name to select the 'probe file
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
  
  
}
