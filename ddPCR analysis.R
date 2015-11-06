# R/SET COMPUTER NAME
input.path <- "D:\\R SCRIPTS\\ddPCR analysis" #work

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
get.breakpoint <- function(x,nClusters=2){ # use kmeans function
  x <- as.numeric(x)
  result <- NULL
  breakpoint <- kmeans(x=x,centers=nClusters)$centers
  if(dim(breakpoint)[1] == 2){result <- mean(breakpoint)}
  if(dim(breakpoint)[1] == 3){result <- c(mean(breakpoint[1:2,1]),mean(breakpoint[2:3,1]))}
  return(result)
}
get.breakpoint.mean <- function(x){ 
  result <- min(x)+((max(x)-min(x))/2)
  return(result)
}
clusters.mean.sd <- function(x,na.rm=TRUE,breakpoint){
    clusters <- c(mean(x[x < breakpoint],na.rm=na.rm), sd(x[x < breakpoint],na.rm=na.rm))
    clusters <- rbind(clusters,c(mean(x[x > breakpoint],na.rm=na.rm), sd(x[x > breakpoint],na.rm=na.rm)))
    clusters <- cbind(clusters,clusters[,2]*3)
    rownames(clusters) <- c("cluster1","cluster2");colnames(clusters) <- c("mean","sd","3*sd")
  return(clusters)
}
create.design <- function(path){
  if(!file.exists(file.path(path,"design.txt"))){
    files <- list.files(path)
    sample.names <- paste("sample",seq(1,length(files),by=1),sep="_")
    sample.type <- c(rep("sample", length(files)-1),"negative")
    output.data <- cbind(Name=sample.names, File=files,Type=sample.type)
    output.file <- file.path(path,"design.txt")
    write.table(x=output.data,file = output.file,quote = FALSE, sep = "\t",row.names = FALSE)
  }
}
# RUN FUNCTIONS
create.folders(path=input.path,c("archive","input.data","output.data","output.plot","scripts","scripts.log"))
path <- set.paths(input.path)
log.file <- paste(format(Sys.time(), "%Y%m%d-%H%M"),"_script_checklist.md",sep="")
make_doc(path=path$scripts,dest = file.path(path$scripts.log,log.file))
# END

# ANALYSIS PER CHANNEL (NOT BOTH AT THE SAME TIME
# RUN THROUGH ALL THE OPTIONS AND CREATE MATRIX FOR EACH CHANNEL
# >10000 DROPLETS = true , INTENSITY>4000= true, ETC... IF ALL OK THEN RUN ANALYSIS, ELSE RESULT IS IN TXT FILE WITH CHECKS.



analysis <- function(){
  # - [x] set input data folder
  data.folder <- file.path(path$input.data,"20151028 EGFR sm_2015-10-28-16-31")
  
  create.design(data.folder)
  # - [x] read input data
  files <- list.files(file.path(path$input.data,data.folder),pattern = ".csv")

  # - [x] create matrix for data check
  nSamples <- length(files)
  results <- matrix(data=NA, nrow=(nSamples*2), ncol=16,byrow = TRUE)
  colnames(results) <- c("FileName","Channel","TotalDroplets","MinAmplitude","MaxAmplitude","DropletsAbove2000Amp","DropletsAbove3000Amp","DropletsAbove4000Amp",
    "Breakpoint","BreakpointKmeans","BreakpointMean","DropletsCluster1_1sd","DropletsCluster2_1sd","DropletsCluster1_3sd","DropletsCluster2_3sd","Overlapping3sdValues")
  counter <- 1
  for(s in 1:length(files)){
    sample.data <- read.table(file.path(path$input.data,data.folder,files[s]),header=TRUE,sep=",")
    result.list <- seq(1,length(files)*2, by=1)
   
    for(c in 1:2){
      if(c == 1){ data <- sample.data[,1]
                  channel <- "Ch1"
                  }
      if(c == 2){ data <- sample.data[,2]
                  channel <- "Ch2" 
                }
      # - [x] find breakpoint between clusters
      data.breakpoint <- get.breakpoint(x=data,nClusters=2)
      data.clusters <- clusters.mean.sd(data,breakpoint=data.breakpoint)
      # - [ ] make a design file for analysis.
      # - [ ] read all the data files and combine.
      # - [ ] find the 2 clusters
      # - [ ] calculate the min/max 3sd of lowest cluster
      # - [ ] find negative samples
      # - [ ] calculate if negative fall between the 3sd of the lowest cluster
      # - [ ] use kmeans as breakpoint for the clusters.
      # - [ ] store breakpoint data for channel 1 & channel 2 for plotting
      
      
      breakpoint.mean <- function(){
        # - [x] find breakpoint between clusters with the mean of Amplitude
        data.breakpoint.mean <- get.breakpoint.mean(x = data)
        data.clusters.mean <- clusters.mean.sd(data,breakpoint=data.breakpoint.mean)
      } # [idea, not yet worked out.]

      # - [x] collect all the results
      results[counter,1] <- files[s]
      results[counter,2] <- channel
      results[counter,3] <- length(data)
      # - [x] what is the max amplitude of the probes
      results[counter,4] <- round(min(data))
      results[counter,5] <- round(max(data))
      # - [x] how much of the probes are located below 2000,3000,4000 intensity?
      results[counter,6] <- sum(data > 2000)
      results[counter,7] <- sum(data > 3000)
      results[counter,8] <- sum(data > 4000)
      results[counter,9] <- data.breakpoint
      results[counter,10] <- round(get.breakpoint(x=data,nClusters=2))
      results[counter,11] <- round(get.breakpoint.mean(x=data))
      # - [x] calculate the amount of probes within 1sd of mean cluater
      results[counter,12] <- sum((data > data.clusters[1,1]-data.clusters[1,2]) & (data < data.clusters[1,1]+data.clusters[1,2])) 
      results[counter,13] <- sum((data > data.clusters[2,1]-data.clusters[2,2]) & (data < data.clusters[2,1]+data.clusters[2,2])) 
      # - [x] calculate the amount of probes within 3sd of mean cluster
      results[counter,14] <- sum((data > data.clusters[1,1]-data.clusters[1,3]) & (data < data.clusters[1,1]+data.clusters[1,3]))
      results[counter,15] <- sum((data > data.clusters[2,1]-data.clusters[2,3]) & (data < data.clusters[2,1]+data.clusters[2,3])) 
      # - [x] calculate if 3sd of mean clusters is overlapping (define the rain cut-offs)
      results[counter,16] <- (data.clusters[2,1]-data.clusters[2,3]) < (data.clusters[1,1]+data.clusters[1,3])
      
      counter <- counter + 1 
    }
    # - [x] create plot for each well
    output.file <- gsub(pattern = ".csv",replacement=".png",x=files[s])
    png(filename = file.path(path$output.plot, output.file),width = 1000,height = 1000)
    plot(sample.data[,1],sample.data[,2], cex=0.5, col="#00000020", ylab="Ch2 Amplitude", xlab="Ch1 Amplitude", pch=19,main=files[s],
         ylim=c(0,11000), xlim=c(0,15000) )
    abline(v=4750, col="red")
    abline(h=4000, col="red")
    dev.off()
  }
  output.file <- file.path(path$output.data,paste(format(Sys.time(), "%Y%m%d"),"_",data.folder,"_ddPCR_analysis.txt",sep=""))
  write.table(file = output.file, x = results,quote = FALSE,sep = "\t",row.names = FALSE)

  
  
  
  # - [ ] [OPTIONAL] Use density lines to find clusters
 	# - [ ] if 2 clusters -> define the rain in amount of droplets
  # - [ ] if 2 clusters -> define the rain in percentage of droplets
 	# - [ ] add column with results to data (0=neg, 1=rain, 2=positive)
 	# - [ ] write data to disk
 	# - [ ] create plot of the data (FAM/HEX)
}
cluster.types <- function(){
  # - [ ] type 1: 2 neg clusters, 1 positive clusters
  # - [ ] type 2: 1 negateive cluster, 1 positive cluster
  # - [ ] type 3: 1 negative cluster
  # - [ ] type 4: 1 negative cluster, 1 positive cluster, noise in between
  # - [ ] type 5: 1 negative cluster, some positive droplets (noise?)
  }
