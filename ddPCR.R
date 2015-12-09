library(dplyr)
library(magrittr)
library(dpcR)
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
  define.color <- function(x,density=NULL)
  {
    ddpcr.colors <- paste(c("#000000","#FF6600","#00CC00","#0033FF"), sep="")
    x <- as.character(x)
    x <- mgsub(pattern = c("1","3","4","2"),replacement = ddpcr.colors,x=x)
    if(class(density) != "NULL")
    {
      x <- paste(x, as.character(density),sep="") 
    }
    
    return(x)
  }
  dropletcount.text <- function(x)
  {
    results <- NULL
    results <- list(clusters=c(cluster.1=sum(x == 1),cluster.2=sum(x == 2),cluster.3=sum(x == 3),cluster.4=sum(x == 4)))
    results  <-paste("Ch1-Ch2-:",results$clusters[1],
                             "   Ch1+Ch2-:",results$clusters[2],
                             "   Ch1+Ch2+:",results$clusters[3],
                             "   Ch1-Ch2+:",results$clusters[4], sep="")
    return(results)
  }
  droplet.count <- function(x,cluster=1)
  {
    result <- sum(x$Cluster %in% cluster)
    return(result)
  }
  get.mean <- function(x,cluster=1,channel=1)
  {
      result <- mean(x[x$Cluster %in% cluster,channel])
      if(as.character(result) == "NaN"){result <- 0}
      return(result)
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
    sub.text <- dropletcount.text(x=x[,3])
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
  cluster.mean.sd <- function(x, cluster = 1, stdev = 1)
  {
    x <- filter(x, Cluster == cluster)
    results <- c(mean(x[,1]), mean(x[,2]))
    results <- rbind(results,c((sd(x[,1])*stdev), (sd(x[,2]))*stdev))
    colnames(results)<- c("Ch1","Ch2")
    rownames(results) <- c("mean",paste(stdev,"_sd",sep=""))
    return(results)
  }
  cluster.range.stdev <- function(x, cluster = 1, stdev = 3)
  { 
    if(cluster == 1)
    {
      result <- colSums(cluster.mean.sd(control.data,cluster = 1,stdev = stdev))
    } else if(cluster == 2)
    {
      result <- cluster.mean.sd(control.data,cluster = 2,stdev = stdev)
      result <- c(Ch1=(result[1,1]-result[2,1]),Ch2=(result[1,2]+result[2,2]))
    } else if(cluster == 3)
    {
      result <- cluster.mean.sd(control.data,cluster = 3,stdev = stdev)
      result <- c(Ch1=(result[1,1]-result[2,1]),Ch2=(result[1,2]-result[2,2]))
    } else if(cluster == 4)
    {
      result <- cluster.mean.sd(control.data,cluster = 4,stdev = stdev)
      result <- c(Ch1=(result[1,1]+result[2,1]),Ch2=(result[1,2]-result[2,2]))
    } else 
    {
      stop("Not the right type of cluster was selected.\n")
    }
  }
  calc.copies <-  function(posCount, count, vDroplet=0.91, volume=1)
  { # concentration in copies / user defined volume
    negCount <- count-posCount
    result <- ((-log(negCount/count)/vDroplet))*1000*volume
    return(result)
  }
  get.plate.wells <- function(prefix=NULL,suffix=NULL)
  {
    result <- rep(NA, 96)
    
    numbers <- c("01","02","03","04","05","06","07","08","09","10","11","12")
    start <- 1
    for(i in 1:8)
    {
      result[start:(start+11)] <- paste(LETTERS[i],numbers, sep="")
      start <- start + 12
    }
    if(class(prefix) != "NULL")
    {
      result <- paste(prefix,result,sep="")
    }
    if(class(suffix) != "NULL")
    {
      result <- paste(result,suffix,sep="")
    }
    return(result)
  }
  get.well <- function(input.file)
  {
    x <- strsplit(as.character(input.file), split = "_")[[1]]
    x <- x[x %in% get.plate.wells()]
    return(x)
  }
  get.statistics <- function(x,sample=NULL,input.file=NULL,target=NULL,breakpoints=NULL, interations=100)
  {
    col.names <- c("Well","Sample","TargetType","Target","Status","CopiesPer1ul","CopiesPer20ulWell","PoissonCopiesPer1ul","ConcentrationPer1ul","Positives",
                   "Negatives","Ch1+Ch2+","Ch1+Ch2-","Ch1-Ch2+","Ch1-Ch2-","AcceptedDroplets","Ratio","FractionalAbundance","Threshold",
                   "MeanAmplitudeofPositives","MeanAmplitudeofNegatives","MeanAmplitudeTotal")
    results <- matrix(NA, nrow = 2, ncol = length(col.names),dimnames = list(NULL,col.names))
    if(class(input.file) != "NULL"){ results[1:2,colnames(results) == "Well"] <- get.well(input.file)}
    if(class(sample) != "NULL"){results[1:2,colnames(results) == "Sample"] <- sample} 
    results[1:2,colnames(results) == "TargetType"] <- c("Ch1","Ch2")
    if(class(target) != "NULL"){ results[1:2,colnames(results) == "Target"] <- sample}
    results[1,colnames(results) == "Positives"] <- droplet.count(x, c(2,3)) #Positives
    results[2,colnames(results) == "Positives"] <- droplet.count(x, c(3,4)) #Positives
    results[1,colnames(results) == "Negatives"] <- droplet.count(x, c(1,4))# Negatives
    results[2,colnames(results) == "Negatives"] <- droplet.count(x, c(1,2)) # Negatives
    results[1:2,colnames(results) == "Ch1+Ch2+"] <- droplet.count(x, 3)
    results[1:2,colnames(results) == "Ch1+Ch2-"] <- droplet.count(x, 2)
    results[1:2,colnames(results) == "Ch1-Ch2+"] <- droplet.count(x, 4)
    results[1:2,colnames(results) == "Ch1-Ch2-"] <- droplet.count(x, 1)
    results[1:2,colnames(results) == "AcceptedDroplets"]<- droplet.count(x, c(1,2,3,4)) #AcceptedDroplets
    if(length(breakpoints) == 2)
    {
      results[1:2,colnames(results) == "Threshold"] <- breakpoints
    }
    results[1,colnames(results) == "CopiesPer1ul"] <- round(calc.copies(posCount = as.numeric(results[1,colnames(results) == "Positives"]),count = as.numeric(results[1,colnames(results) == "AcceptedDroplets"])), digits = 1)
    results[2,colnames(results) == "CopiesPer1ul"] <- round(calc.copies(posCount = as.numeric(results[2,colnames(results) == "Positives"]),count = as.numeric(results[2,colnames(results) == "AcceptedDroplets"])), digits = 1)
    results[1,colnames(results) == "PoissonCopiesPer1ul"] <- round(calc.copies(posCount = poisson.correction(count = as.numeric(results[1,colnames(results) == "AcceptedDroplets"]),
                                                                                                                posCount = as.numeric(results[1,colnames(results) == "Positives"]), iterations = interations),
                                                                                                                count = as.numeric(results[1,colnames(results) == "AcceptedDroplets"])), digits = 1)
    results[2,colnames(results) == "PoissonCopiesPer1ul"] <- round(calc.copies(posCount = poisson.correction(count = as.numeric(results[2,colnames(results) == "AcceptedDroplets"]),
                                                                                                                posCount = as.numeric(results[2,colnames(results) == "Positives"]), iterations = interations),
                                                                                                                count = as.numeric(results[2,colnames(results) == "AcceptedDroplets"])), digits = 1)
    results[1:2,colnames(results) == "CopiesPer20ulWell"] <- as.numeric(results[1:2,colnames(results) == "CopiesPer1ul"])*20
    results[1:2,colnames(results) == "ConcentrationPer1ul"] <- as.numeric(results[1:2,colnames(results) == "CopiesPer1ul"])/15.152
    results[1:2,colnames(results) == "Ratio"] <- as.numeric(results[1,colnames(results) == "CopiesPer1ul"]) / as.numeric(results[2,colnames(results) == "CopiesPer1ul"])
    results[2,colnames(results) == "Ratio"] <- as.numeric(results[2,colnames(results) == "CopiesPer1ul"]) / as.numeric(results[1,colnames(results) == "CopiesPer1ul"])
    results[1,colnames(results) == "FractionalAbundance"] <- as.numeric(results[1,colnames(results) == "CopiesPer1ul"]) / (as.numeric(results[1,colnames(results) == "CopiesPer1ul"])+as.numeric(results[2,colnames(results) == "CopiesPer1ul"]))*100
    results[2,colnames(results) == "FractionalAbundance"] <- as.numeric(results[2,colnames(results) == "CopiesPer1ul"]) / (as.numeric(results[1,colnames(results) == "CopiesPer1ul"])+as.numeric(results[2,colnames(results) == "CopiesPer1ul"]))*100
    results[1,colnames(results) == "MeanAmplitudeofPositives"] <- round(get.mean(x,cluster=c(2,3),1), digits = 2)
    results[2,colnames(results) == "MeanAmplitudeofPositives"] <- round(get.mean(x,cluster=c(4,3),2), digits = 2)
    results[1,colnames(results) == "MeanAmplitudeofNegatives"] <- round(get.mean(x,cluster=c(1,4),1), digits = 2)
    results[2,colnames(results) == "MeanAmplitudeofNegatives"] <- round(get.mean(x,cluster=c(1,2),2), digits = 2)
    results[1,colnames(results) == "MeanAmplitudeTotal"] <- round(get.mean(x,cluster=c(1,2,3,4),1), digits = 2)
    results[2,colnames(results) == "MeanAmplitudeTotal"] <- round(get.mean(x,cluster=c(1,2,3,4),2), digits = 2)
    results[1:2,colnames(results) == "Status"] <- "OK"
    if(as.numeric(results[1:2,colnames(results) == "AcceptedDroplets"]) < 10000){results[1:2,colnames(results) == "Status"] <- "CHECK"}
    if(as.numeric(results[1,colnames(results) == "CopiesPer1ul"]) < 10000){results[1:2,colnames(results) == "Status"] <- "CHECK"}
    if(as.numeric(results[2,colnames(results) == "CopiesPer1ul"]) < 10000){results[1:2,colnames(results) == "Status"] <- "CHECK"}
    return(results)
  }
  poisson.correction <- function(posCount, count, iterations=100)
  {
    if(posCount == 0)
    {
      return(0)
    } else 
      {
        additional.droplets <- rep(NA,iterations)
        for(i in 1:iterations)
        {
          result <- rep(NA,count)
          droplets <- 1:count
          result[1:posCount] <- sample(x = droplets,size = posCount,replace = TRUE)
          counter <- posCount+1
          for(z in (posCount+1):count)
          {
            if((length(unique(result))-1) >= posCount)
            {
              break
            }
            result[z] <- sample(x = droplets,size = 1,replace = TRUE)
            z <- z+1
          }
          additional.droplets[i] <- z-posCount
        }
        additional.droplets <- mean(c((posCount+min(additional.droplets)),posCount+max(additional.droplets)))
        return(additional.droplets)
      }
  }
  convert.copies.to.ng <- function(x)
  {
    x <- x / 15.152
    return(x)
  }
  exactPoiCI <- function (x, conf.level=0.95) 
  {
  poisson.test(x, conf.level = conf.level)$conf.int[1:2]
  }
  read.csv <- function(file)
  {
    temp <- read.table(file = file,header = FALSE,sep = "\t",fill=TRUE)
    temp <- as.matrix(temp)
    nr.lines <- dim(temp)[1]
    nr.cols <- max(apply(X = temp,MARGIN = 1,FUN=countCharOccurrences))+1
    test.matrix <- matrix("NA",ncol=nr.cols,nrow=nr.lines)
    for (d in 1:(nr.lines)){
      x <-  unlist(strsplit(temp[d,1],split = ","))
      columns <- length(x)
      test.matrix[d,1:columns] <- x[1:columns]
    }
    return(test.matrix)
  }
  get.experiments <- function(path)
  {
    file <- list.files(path, pattern = "Error.log",full.names = TRUE)
    file <- gsub(pattern = "Error.log", replacement = "",x = file)
    experiment.name <- basename(file)
    file <- paste(file,".csv", sep="")
    if(file.exists(file) == TRUE)
    {
      data <-  read.table(file = file,header = TRUE, sep=",",row.names = NULL)
    } else { break 
    }
    targets <- unique(gsub(pattern = " wt", replacement = "", x=data$Target))
    grep(data$Target,pattern = targets[i])
    
    result <- list()
    for(i in 1:length(targets))
    {
      result[[i]] <- data[grep(data$Target,pattern = targets[i]),]
    }
    names(result) <- targets
    return(result)
  }
  
  # end, HF van Essen 2015