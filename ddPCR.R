library(dplyr)
library(magrittr)
# FUNCTIONS
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
  } # currently not used
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
  } # currently not used
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
  {
    x <- as.numeric(x)
    result <- NULL
    breakpoint <- kmeans(x=x,centers=nClusters)$centers
    if(dim(breakpoint)[1] == 2){result <- mean(breakpoint)}
    return(result)
  }
  get.breakpoint.ranges <- function(x)
  {
    x <- as.numeric(x)
    result <- NULL
    result <-  (max(x) - min(x))/2
    return(result)
  }
  get.breakpoint.hist <- function(x)
  {
    x <- as.numeric(x)
    hist.data <- rbind(hist(x, breaks=15, plot=FALSE)$mids, hist(x, breaks=15, plot=FALSE)$counts)
    hist.data <- hist.data[,-c(1:2,14:16)]
    result <- mean(hist.data[1,hist.data[2,] == min(hist.data[2,])])
    return(result)
  }
  get.ddpcr.breakpoints.kmeans <- function(x)
  {
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
  plot.cutoffs <- function(x,col="black")
  {
    abline(h=x[1], col=col)
    abline(v=x[2], col=col)
  }
  mean.cluster <- function(x,cluster=1, zero=TRUE)
  {
    if(zero == TRUE)
    {
      x <- c(0,0);
      x <- rbind(x,c(0,0))
      colnames(x)<- c("Ch1","Ch2")
      rownames(x) <- c("mean","sd")
    } else 
    {
      x %<>% 
        data.frame(.) %>%
        filter(.,Cluster == cluster) %>%
        cluster.mean.sd(.)
    }
    return(x)
  }
  cluster.mean.sd <- function(x)
  {
    cluster <- c(mean(x[,1]), mean(x[,2]))
    cluster <- rbind(cluster,c(sd(x[,1]), sd(x[,2])))
    colnames(cluster)<- c("Ch1","Ch2")
    rownames(cluster) <- c("mean","sd")
    return(cluster)
  }
  mean.4sd.cutoff <- function(x)
  { # mean of the clusters + 4 times the sd of the cluster
    cluster.1 <- mean.cluster(x, cluster=1)
    if(sum(x$Cluster == 2) == 0)
      {
        cluster.2 <- mean.cluster(x, cluster=2, zero=TRUE)
      } else
      {
        cluster.2 <- mean.cluster(x, cluster=2)
      }
    if(sum(x$Cluster == 4) == 0)
    {
      cluster.4 <- mean.cluster(x, cluster=4, zero=TRUE)
    } else
    {
      cluster.4 <- mean.cluster(x, cluster=4)
    }
    channel.1 <- max(cluster.1[1,1]+(cluster.1[2,1]*4),cluster.4[1,1]+(cluster.4[2,1]*4))
    channel.2 <- max(cluster.1[1,2]+(cluster.1[2,2]*4),cluster.2[1,2]+(cluster.2[2,2]*4))
    result <- c(Ch1=channel.1,Ch2=channel.2)
    return(result)  
  }
# end, HF van Essen 2015