

#The concentration is then calculated using the following formula:
 c = -ln ((Nneg /N)Vdroplet)
  ###

 
 
### use histo to determine the two clusters
sample <- hist(train,breaks=2)
breaks <- train
### determine the correct means of each group based on the breaks

############### define my rain
select.file.path <- function(){
   location <- file.choose()
   inputFile <- basename(location)
   pad <- unlist(strsplit(location, inputFile))
   setwd(pad)
   data <- c(inputFile,pad)
   return(data)
 }
parse.csv <- function(input.file){
   nr.lines       <- length(count.fields(input.file,blank.lines.skip = FALSE))
   test.matrix <- matrix("NA",ncol=15,nrow=nr.lines)
   for (d in 1:(nr.lines)){
     temp.data <-  read.csv(input.file,header=FALSE,sep=",", fill=TRUE,skip=d-1,blank.lines.skip = TRUE,nrows=1,na.strings = "NA") 
     for(i in 1:length(temp.data)){test.matrix[d,i] <- as.character(temp.data[1,i])}
   }
   return(test.matrix)
 }
clusters.mean.sd <- function(data,na=TRUE,breakpoint){
  neg.mean <- mean(x[x < breakpoint],na.rm=na)
  neg.sd <- sd(x[x < breakpoint],na.rm=na)
  pos.mean <- mean(x[x > breakpoint],na.rm=na)
  pos.sd <- sd(x[x > breakpoint],na.rm=na)
  return(list(neg.mean=neg.mean,neg.sd=neg.sd,pos.mean=pos.mean,pos.sd=pos.sd))
}
get.clusters <- function(x){
  x <- as.numeric(x)
  hist.data <- hist(x,breaks=2)
  if (length(hist.data$mids) == 2){
    result <- clusters.mean.sd(data=x,na=TRUE,breakpoint=hist.data$breaks[2])
  }
  if (length(temp$mids) == 3){
    result <- clusters.mean.sd(data=x,na=TRUE,breakpoint=hist.data$mids[2])
  }
  return(result)
}

data.mean.sd <- get.clusters(x=x)

plot(data[,1])
 
 ######## 2 sets with rain
file <- select.file.path()
data <- read.table(file[1],",",header=TRUE)
data.mean.sd <- get.clusters(x=data[,1])
rain <- x > (data.mean.sd$neg.mean + (3 * data.mean.sd$neg.sd)) & x < (data.mean.sd$pos.mean - (3 * data.mean.sd$pos.sd))

color <- rep("black",length(rain));color[rain] <- "red";color <- as.character(color)
plot(x,cex=0.3,col=color)
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
