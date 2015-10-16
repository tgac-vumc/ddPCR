############################################
#
#             2 CHANNEL, 2 CLUSTERS
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
    clusters <- cbind(clusters,clusters[,2]*2,clusters[,2]*3)
    rownames(clusters) <- c("cluster1","cluster2");colnames(clusters) <- c("mean","sd","2*sd","3*sd")
  }
  if(length(breakpoints) == 2){
    breakpoints <- sort(breakpoints)
    clusters <- c(mean(x[x < breakpoints[1]],na.rm=na), sd(x[x < breakpoints[1]],na.rm=na))
    clusters <- rbind(clusters,c(mean(x[x > breakpoints[1] & x < breakpoints[2]],na.rm=na), sd(x[x > breakpoints[1] & x < breakpoints[2]],na.rm=na)))
    clusters <- rbind(clusters,c(mean(x[x > breakpoints[2]],na.rm=na), sd(x[x > breakpoints[2]],na.rm=na)))
    clusters <- cbind(clusters,clusters[,2]*2,clusters[,2]*3)
    rownames(clusters) <- c("cluster1","cluster2","cluster3");colnames(clusters) <- c("mean","sd","2*sd","3*sd")
  }
  return(clusters)
}
calculate.cutoff <- function(x,sd.value){
  result <- matrix(data = "NA",nrow = 1, ncol = 2)
  result[1,1] <- x-sd.value
  result[1,2] <- x+sd.value
  colnames(result) <- c("cutoff_low","cutoff_high")
  return(result)
}
######
if(as.character(Sys.info()["nodename"]) == "PA"){
  setwd("C:\\Documents and Settings\\h.vanessen\\Mijn documenten\\Dropbox\\r scripts\\define_the_rain\\data\\other\\")
}
if(as.character(Sys.info()["nodename"]) == "MacBook-Air-van-Dirk.local"){
  setwd("/Users/dirkvanessen/.dropbox-two/Dropbox/r scripts/define_the_rain/DATA/other/")
}

files <- list.files(pattern = ".csv")
result.ch1 <- matrix("NA", ncol=6, nrow=length(files))
  colnames(result.ch1) <- c("file","channel","pos","neg","rain","warnings")
result.ch2 <- matrix("NA", ncol=6, nrow=length(files))
  colnames(result.ch2) <- c("file","channel","pos","neg","rain","warnings")
for(i in 1:length(files)){
  data <- read.table(files[i],",",header=TRUE)
  nClusters=2
  #plot(data,cex=1,pch=20)
  for(c in 1:2){
    if(c == 1){data <- data[order(data[,2]),]
      }else{data <- data[order(data[,1]),]
    }
    dye <- colnames(data[c])
    breakpoints <- get.breakpoints(x=data[,c],nClusters=nClusters)
    clusters <- clusters.mean.sd(data[,c],breakpoints=breakpoints)
    if(c == 1){clusters <- calculate.cutoff(x = clusters,sd.values=clusters[2,4])
              }else{clusters <- calculate.cutoff(x = clusters,sd.values=clusters[1,4])
    }
### I AM HERE
    colors <- rep("blue",length(data[,c]))
    colors[data[,c] < clusters[1,5]] <- "red"; 
    colors[data[,c] > clusters[2,4]] <- "green"
    results <- data.frame(sum(colors == "green"),sum(colors == "red"),sum(colors == "blue"));
    colnames(results) <- c("pos","neg","rain")
      if(c == 1){ result.ch1[i,1] <- files[i]
                  result.ch1[i,2] <- dye
                  result.ch1[i,3:5] <- as.character(results)
      }
      if(c == 2){ result.ch2[i,1] <- files[i]
                  result.ch2[i,2] <- dye
                  result.ch2[i,3:5] <- as.character(results)
      }
    sub.text <- paste("positive=",results[1],"; negative=",results[2],"; rain=",results[3])
    if (!file.exists('profiles'))
      dir.create('profiles')
    output.file <- paste(gsub("Amplitude.csv","",files[i]),dye,".png",sep="")
    png(filename = file.path('profiles',output.file))
    plot(data[,c],cex=1,pch=20,col=colors, main=dye,ylab="Intensity",xlab="nProbes")
    mtext(sub.text,side = 3,adj = 1)
    dev.off()
  }
}
save.date.time <- format(Sys.time(), "%Y%m%d-%H%M")
output.file <- paste(save.date.time,"_ddPCR_results.txt",sep="")
output.data <- rbind(result.ch1,result.ch2)
write.table(x = output.data, file = output.file, row.names = FALSE, quote = FALSE,sep="\t")
