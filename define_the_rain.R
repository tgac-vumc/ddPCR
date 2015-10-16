### DEFINE THE RAIN
archive <- function(){
  #### The concentration is then calculated using the following formula:
 #   c = -ln ((Nneg /N)Vdroplet)
  ### use histo to determine the two clusters
  sample <- hist(train,breaks=2)
  breaks <- train
  ### determine the correct means of each group based on the breaks
  get.clusters <- function(x){ ###old functionS
    x <- as.numeric(x)
    result <- NULL
    hist.data <- hist(x,breaks=2)
    if (length(hist.data$mids) == 2){
      result <- clusters.mean.sd(x=x,na=TRUE,breakpoint=hist.data$breaks[2])
    }
    if (length(hist.data$mids) == 3){
      result <- clusters.mean.sd(x=x,na=TRUE,breakpoint=hist.data$mids[2])
    }
    if(class(result) == "NULL"){cat("Error in 'get.clusters()': breaks not properly detected.\n")
    }else{return(result)}
  }
  #### from get.clusters function when more clusters are needed
  if (length(hist.data$mids) == 3){
    result <- clusters.mean.sd(x=x,na=TRUE,breakpoint=hist.data$mids[2])
  }
}
############### define my rain FUNCTIONS
select.file.path <- function(){
  location <- file.choose()
  inputFile <- basename(location)
  pad <- unlist(strsplit(location, inputFile))
  setwd(pad)
  data <- c(inputFile,pad)
  return(data)
}
clusters.mean.sd <- function(x,na=TRUE,breakpoint){
  neg.mean <- mean(x[x < breakpoint],na.rm=na)
  neg.sd <- sd(x[x < breakpoint],na.rm=na)
  pos.mean <- mean(x[x > breakpoint],na.rm=na)
  pos.sd <- sd(x[x > breakpoint],na.rm=na)
  return(list(neg.mean=neg.mean,neg.sd=neg.sd,pos.mean=pos.mean,pos.sd=pos.sd))
}
get.clusters <- function(x,nClusters=2){ # use kmeans function
  x <- as.numeric(x)
  result <- NULL
  breakpoint <- mean(kmeans(x=x,centers=nClusters)$centers)
  if (nClusters == 2){
    result <- clusters.mean.sd(x=x,na=TRUE,breakpoint=breakpoint)
  }
  if(class(result) == "NULL"){cat("Error in 'get.clusters()': breaks not properly detected.\n")
  }else{return(result)}
}
define.rain.old <- function(x,data.mean.sd){
  rain <- x > (data.mean.sd$neg.mean + (3 * data.mean.sd$neg.sd)) & x < (data.mean.sd$pos.mean - (3 * data.mean.sd$pos.sd))
  return(rain)
}
define.rain.new <- function(x,data.mean.sd){
  
  rain <- x > (data.mean.sd$neg.mean + (3 * data.mean.sd$neg.sd)) & x < (data.mean.sd$pos.mean - (3 * data.mean.sd$pos.sd))
  return(rain)
}
##### END OF FUNCTIONS

   ######## 2 channels, 4 clusters
if(as.character(Sys.info()["nodename"]) == "PA"){
  setwd("C:\\Documents and Settings\\h.vanessen\\Mijn documenten\\Dropbox\\r scripts\\define_the_rain\\data\\other\\")
}
if(as.character(Sys.info()["nodename"]) == "MacBook-Air-van-Dirk.local"){
  setwd("/Users/dirkvanessen/.dropbox-two/Dropbox/r scripts/define_the_rain/DATA/other/")
}

files <- list.files(pattern = ".csv")
for(i in 1:length(files)){
  data <- read.table(files[i],",",header=TRUE)
  data.mean.sd <- get.clusters(x=data[,1])
  rain1 <- define.rain(x=data[,1],data.mean.sd=data.mean.sd) 
  data.mean.sd <- get.clusters(x=data[,2])
  rain2 <- define.rain(x=data[,2],data.mean.sd=data.mean.sd) 
  
  color <- rep("black",length(rain));
  color[rain1] <- "red";
  color[rain2] <- "red";
  color <- as.character(color)
  
  ouput.file <- gsub(pattern = ".csv",replacement = ".png",files[i])
  png(filename = ouput.file)
  plot(data,cex=1,col=color,pch=20)
  dev.off()
};i=1

############################################
#
#             2 CLUSTERS RESTART
#
############################################
get.breakpoints <- function(x,nClusters=2){ # use kmeans function
  x <- as.numeric(x)
  result <- NULL
  breakpoint <- kmeans(x=x,centers=nClusters)$centers
  if(dim(breakpoint)[1] == 2){result <- mean(breakpoint)}
  if(dim(breakpoint)[1] == 3){result <- c(mean(breakpoint[1:2,1]),mean(breakpoint[2:3,1]))}
  return(result)
}
clusters.mean.sd <- function(x,na=TRUE,breakpoints){
  if(length(breakpoints) == 1){
   clusters <- c(mean(x[x < breakpoints],na.rm=na), sd(x[x < breakpoints],na.rm=na))
   clusters <- rbind(clusters,c(mean(x[x > breakpoints],na.rm=na), sd(x[x > breakpoints],na.rm=na)))
   clusters <- cbind(clusters,clusters[,2]*3)
   rownames(clusters) <- c("cluster1","cluster2");colnames(clusters) <- c("mean","sd","3*sd")
  }
  if(length(breakpoints) == 2){
    breakpoints <- sort(breakpoints)
    clusters <- c(mean(x[x < breakpoints[1]],na.rm=na), sd(x[x < breakpoints[1]],na.rm=na))
    clusters <- rbind(clusters,c(mean(x[x > breakpoints[1] & x < breakpoints[2]],na.rm=na), sd(x[x > breakpoints[1] & x < breakpoints[2]],na.rm=na)))
    clusters <- rbind(clusters,c(mean(x[x > breakpoints[2]],na.rm=na), sd(x[x > breakpoints[2]],na.rm=na)))
    clusters <- cbind(clusters,clusters[,2]*3)
    rownames(clusters) <- c("cluster1","cluster2","cluster3");colnames(clusters) <- c("mean","sd","3*sd")
  }
  return(clusters)
}
calculate.cutoff <- function(x){
  x <- cbind(x,(x[,1]-x[,3]))
  x <- cbind(x,(x[,1]+x[,3]))
  colnames(x)[4:5] <- c("cutoff_low","cutoff_high")
  return(x)
}
######
 if(as.character(Sys.info()["nodename"]) == "PA"){
    setwd("C:\\Documents and Settings\\h.vanessen\\Mijn documenten\\Dropbox\\r scripts\\define_the_rain\\data\\other\\")
  }
 if(as.character(Sys.info()["nodename"]) == "MacBook-Air-van-Dirk.local"){
    setwd("/Users/dirkvanessen/.dropbox-two/Dropbox/r scripts/define_the_rain/DATA/other/")
  }
  
  files <- list.files(pattern = ".csv")
  for(i in 1:length(files)){
  };i=1
    data <- read.table(files[i],",",header=TRUE)
    nClusters=2
    plot(data,cex=1,pch=20)
    for(c in 1:2){
      breakpoints <- get.breakpoints(x=data[,c],nClusters=nClusters)
      clusters <- clusters.mean.sd(data[,c],breakpoints=breakpoints)
      clusters <- calculate.cutoff(x = clusters)
      if(c == 1){abline(v=clusters[,4:5],col="red")}
      if(c == 2){abline(h=clusters[,4:5],col="red")}
      print(clusters)
    }
    
  
 # 1. selection 2 positive
#  2. selection rain
#  3. selection 2 negative

############################################
#
#             2 CLUSTERS END
#
############################################


### TO DO LIST
#manually set the amount of clusters (2 or 3)
#get mean and sd for each cluster
#for 2 clusters 
#  determine high and low cluster
#  high == mean - 3*sd
#  low == mean + 3*sd
  
#adjust rain output data for each channel
#  - neg / rain cut-off
#  - pos / rain cut-off
#create color vector based on rain.output
#  - red  / double neg
#  - green / double pos
#  - blue / one channel pos
#  - black / rain

#create output list$double.pos, list$pos1, list$pos2, list$neg,
####
three.clusters <- function(){
  ### 3 sets with neg and sample and adapters
  x <- data[,1]
  hist.data <- hist(x,breaks=3)
  color <- rep("black",length(x));
  color[x < hist.data$breaks[3]] <- "red"
  color[x > hist.data$breaks[5]] <- "blue"
  color <- as.character(color)
  plot(x,cex=0.3,col=color)
  # test sets
  x <- c(rnorm(mean = 5,n = 500),rnorm(n=300,mean=30)) # 2 clusters
  x <- c(rnorm(mean = 5,n = 5000),rnorm(n=3000,mean=10),rnorm(n=100,mean = 15)) # 2 clusters
  ###
}
plot.all.csv <- function(){
  #### plot all csv files
  files <- list.files(select.file.path()[2])
  
  for(i in 1:length(files)){
    data <- read.table(files[i],header=TRUE,sep=",")
    output.file <- paste(gsub(pattern = ".csv",replacement = "",files[i]),".png",sep="")
    png(output.file)
    plot(data[,1],data[,2])
    dev.off()
  }
}


 

