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
    archive <- file.path(path,"archive")
    input.data <-  file.path(path,"input.data")
    output.data <-  file.path(path,"output.data")
    scripts <-  file.path(path,"scripts")
    scripts.log <- file.path(path,"scripts.log")
    # - [x] set paths
    paths <- list(archive=archive,
                  input.data=input.data,
                  output.data=output.data,
                  scripts=scripts,
                  scripts.log=scripts.log)
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
clusters.mean.sd <- function(x,na.rm=TRUE,breakpoint){
    clusters <- c(mean(x[x < breakpoint],na.rm=na.rm), sd(x[x < breakpoint],na.rm=na.rm))
    clusters <- rbind(clusters,c(mean(x[x > breakpoint],na.rm=na.rm), sd(x[x > breakpoint],na.rm=na.rm)))
    clusters <- cbind(clusters,clusters[,2]*3)
    rownames(clusters) <- c("cluster1","cluster2");colnames(clusters) <- c("mean","sd","3*sd")
  return(clusters)
}
# RUN FUNCTIONS
create.folders(path=input.path,c("archive","input.data","output.data","scripts","scripts.log"))
path <- set.paths(input.path)
log.file <- paste(format(Sys.time(), "%Y%m%d-%H%M"),"_script_checklist.md",sep="")
make_doc(path=path$scripts,dest = file.path(path$scripts.log,log.file))
# END

analysis <- function(){
  # - [x] set input data folder
  data.folder <- "20150226 nano quantification"
  # - [x] read input data
  files <- list.files(file.path(path$input.data,data.folder),pattern = ".csv")
  i=11
  data <- read.table(file.path(path$input.data,data.folder,files[i]),header=TRUE,sep=",")
  plot(data[,1],data[,2], cex=0.5, col="00000020", ylab="Channel 1", xlab="Channel 2")
 	# - [x] check amount of droplets
  if(dim(data)[1] < 10000){
    cat("Number of droplets for well", files[i], "is not enough for automated analysis (< 10.000). Sample will be skipped.\n")
    next
  }else{cat("A total of", dim(data)[1], "droplets were found for analysis.\n")}
	# - [x] check the intensity of the droplets
  if(max(data[,1]) < 4000){
    cat("Intensity for channel 1 for well", files[i], "is lower then expected. Please check your sample.\n")
  } else {cat("Maximum intensity for channel 1 for well ", files[i], " is ",max(data[,1]),".\n", sep="")
    }
  if(max(data[,1]) < 4000){
    cat("Intensity for channel 2 for well", files[i], "is lower then expected. Please check your sample.\n")
  } else {cat("Maximum intensity for channel 2 for well ", files[i], " is ",max(data[,1]),".\n",sep="")
    }
	# - [x] find clusters (1, 2, 3, 4)
  breakpoints.ch1 <- get.breakpoints(x=data[,1],nClusters=2)
  breakpoints.ch2 <- get.breakpoints(x=data[,2],nClusters=2)
  clusters.ch1 <- clusters.mean.sd(data[,1],breakpoints=breakpoints.ch1)
  clusters.ch2 <- clusters.mean.sd(data[,2],breakpoints=breakpoints.ch2)
  # - [ ] [OPTIONAL] Use density lines to find clusters
  
 	# - [ ] if 2 clusters -> define the rain (FAM/HEX)
 	# - [ ] add column with results to data (0=neg, 1=rain, 2=positive)
 	# - [ ] write data to disk
 	# - [ ] create plot of the data (FAM/HEX)
} 