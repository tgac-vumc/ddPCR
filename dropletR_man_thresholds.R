##########################################################################################
#   TITLE : drops.pl : Data processing for droplet digital PCR data (ddPCR v1.01)        #
##########################################################################################
##########################################################################################
#									     												 #
#	This R script is designed to process data from the Bio-Rad DX100 platform		     #
#	For details, see the perl script "drops.pl"											 #
#																						 #
#	For most purposes, the concise report file is sufficient to get the job done.		 #
##########################################################################################
#                                   Open source license                                  #
##########################################################################################
# To attribute this work, you must cite the name of the original source, the name of     #
# the author (Chrissy h. Roberts) and their contact details                              #
# (chrissyhroberts@yahoo.co.uk). This work is licensed under the Creative Commons        #
# Attribution-ShareAlike 3.0 Unported License.                                           #
# To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/3.0/   #
# or send a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,     #
# California, 94041, USA.                                                                #
##########################################################################################
#
##########################################################################################
#START OF SCRIPT
##########################################################################################


delete.all <- function()
     rm(list=ls(pos=.GlobalEnv), pos=.GlobalEnv)
  delete.all()
ls()
cwd <-getwd()
cwd
cat(cwd,"/out/",sep="")

outwd <-paste(cwd,"/out/",sep="")

path <- cwd

outpath <- outwd

#READ IN PLATE
countCharOccurrences <- function(x,header=FALSE) {
  nchar.start <- nchar(x)
  nchar.end <- nchar(gsub(x,pattern = ",",replacement = ""))
  return (nchar.start - nchar.end)
}
read.csv <- function(file,header=FALSE){
  temp <- read.table(file = file,sep = "\t",fill=TRUE)
  temp <- as.matrix(temp)
  nr.lines <- dim(temp)[1]
  nr.cols <- max(apply(X = temp,MARGIN = 1,FUN=countCharOccurrences))+1
  test.matrix <- matrix("NA",ncol=nr.cols,nrow=nr.lines)
  for (d in 1:(nr.lines)){
    x <-  unlist(strsplit(temp[d,1],split = ","))
    columns <- length(x)
    test.matrix[d,1:columns] <- x[1:columns]
  }
  if(header == TRUE){
    names <- test.matrix[1,]
    test.matrix <- test.matrix[2:dim(test.matrix)[1],]
    colnames(test.matrix) <- names
  }
  return(test.matrix)
}

select.file.path <- function(){
  location <- file.choose()
  inputFile <- basename(location)
  pad <- unlist(strsplit(location, inputFile))
  setwd(pad)
  data <- c(inputFile,pad)
  return(data)
}
path <- "/Users/dirkvanessen/.dropbox-two/Dropbox/ddPCR sample calculation/"
setwd(path)
file <- "Test Illumina Quantification.csv"
PLATE <- read.csv(file=file,header=TRUE)


#PLOT TOTAL DATA FAM VS VIC

pdf("platesummary.pdf")
plot(PLATE$FAM,PLATE$VIC,xlab="PLASMID",ylab="HURNASE", xlim=c(0,max(PLATE$FAM)+100),ylim=c(0,max(PLATE$VIC)+100),pch = 46)
dev.off()


samples<-table(PLATE$SAMPLE)
samples<-as.data.frame(samples)
samples$chlamydia_test<-NA
samples$FAILS<-NA
samples$concentration_fam_per_swab<-NA
samples$concentration_fam_per_swab_low<-NA
samples$concentration_fam_per_swab_high<-NA
samples$concentration_vic_per_swab<-NA
samples$concentration_vic_per_swab_low<-NA
samples$concentration_vic_per_swab_high<-NA

rownumbers <- c(1:length(samples$Var1))


#sample by sample thresholds and gates

rownumbers <- c(1:length(samples$Var1))

for (i in rownumbers) {

#get well name
coordinate <- samples$Var1[i]

#subset based on well name
currentsubset <- subset(PLATE,PLATE$SAMPLE==coordinate)
x11()
plot(PLATE$FAM[PLATE$SAMPLE==samples$Var1[i]],PLATE$VIC[PLATE$SAMPLE==samples$Var1[i]],xlab="PLASMID",ylab="HURNASE", xlim=c(0,max(PLATE$FAM)+100),ylim=c(0,max(PLATE$VIC)+100),pch = 46,main=samples$Var1[i])
thresholds<-locator(n=1)
dev.off()

samples$famthreshold[i]<-thresholds$x[1]
samples$victhreshold[i]<-thresholds$y[1]

#count fam pos only
samples$famposvicneg[i]<-((length(which(currentsubset$FAM > samples$famthreshold[i] & currentsubset$VIC < samples$victhreshold[i]) )))
if (is.na(samples$famposvicneg[i])){samples$famposvicneg[i]<-0}


#count vic pos only
samples$famnegvicpos[i]<-((length(which(currentsubset$FAM < samples$famthreshold[i] & currentsubset$VIC > samples$victhreshold[i]) )))
if (is.na(samples$famnegvicpos[i])){samples$famnegvicpos[i]<-0}

#count double pos only
samples$famposvicpos[i]<-((length(which(currentsubset$FAM > samples$famthreshold[i] & currentsubset$VIC > samples$victhreshold[i]) )))
if (is.na(samples$famposvicpos[i])){samples$famposvicpos[i]<-0}

#count double neg only
samples$famnegvicneg[i]<-((length(which(currentsubset$FAM < samples$famthreshold[i] & currentsubset$VIC < samples$victhreshold[i]) )))
if (is.na(samples$famnegvicneg[i])){samples$famnegvicneg[i]<-0}

pdf(paste(outpath,samples$Var1[i], ".pdf", sep=""))
print (samples$Var1[i])
plot(PLATE$FAM[PLATE$SAMPLE==samples$Var1[i]],PLATE$VIC[PLATE$SAMPLE==samples$Var1[i]],xlab="PLASMID",ylab="HURNASE", xlim=c(0,max(PLATE$FAM)+100),ylim=c(0,max(PLATE$VIC)+100),pch = 46,main=samples$Var1[i])
abline(v=samples$famthreshold[i], col="PINK", lty=3)
abline(h=samples$victhreshold[i], col="PINK", lty=3)
dev.off()



#CALCULATE AVERAGE POS FI AND AVERAGE NEG FI


currentsubset_2<-subset(currentsubset,(currentsubset$FAM > samples$famthreshold[i]))
samples$ave_fi_fam_positives[i]<-mean(currentsubset_2$FAM)


currentsubset_2<-subset(currentsubset,(currentsubset$FAM < samples$famthreshold[i]))
samples$ave_fi_fam_negatives[i]<-mean(currentsubset_2$FAM)


currentsubset_2<-subset(currentsubset,(currentsubset$VIC > samples$victhreshold[i]))
samples$ave_fi_vic_positives[i]<-mean(currentsubset_2$VIC)



currentsubset_2<-subset(currentsubset,(currentsubset$VIC < samples$victhreshold[i]))
samples$ave_fi_vic_negatives[i]<-mean(currentsubset_2$VIC)





#remove currentsubset for tidiness

rm(currentsubset,currentsubset_2)
}

 # quality assessment of FI values from positives and negatives

platewide_ave_fam_pos_FI<-mean(samples$ave_fi_fam_positives,na.rm=T)
platewide_ave_fam_neg_FI<-mean(samples$ave_fi_fam_negatives,na.rm=T)
platewide_ave_vic_pos_FI<-mean(samples$ave_fi_vic_positives,na.rm=T)
platewide_ave_vic_neg_FI<-mean(samples$ave_fi_vic_negatives,na.rm=T)
sdplatewide_ave_fam_pos_FI<-sd(samples$ave_fi_fam_positives,na.rm=T)
sdplatewide_ave_fam_neg_FI<-sd(samples$ave_fi_fam_negatives,na.rm=T)
sdplatewide_ave_vic_pos_FI<-sd(samples$ave_fi_vic_positives,na.rm=T)
sdplatewide_ave_vic_neg_FI<-sd(samples$ave_fi_vic_negatives,na.rm=T)
if (is.na(sdplatewide_ave_fam_neg_FI)){sdplatewide_ave_fam_neg_FI<-platewide_ave_fam_neg_FI}
if (is.na(sdplatewide_ave_fam_pos_FI)){sdplatewide_ave_fam_pos_FI<-platewide_ave_fam_pos_FI}
if (is.na(sdplatewide_ave_vic_neg_FI)){sdplatewide_ave_vic_neg_FI<-platewide_ave_vic_neg_FI}
if (is.na(sdplatewide_ave_vic_pos_FI)){sdplatewide_ave_vic_pos_FI<-platewide_ave_vic_pos_FI}

#plot vic average FI for positives and negatives with platewide 99%CI
pdf(file.path(outpath,"000_average_vic_FI.pdf"))

plot(samples$ave_fi_vic_positives,pch=46,col="red",cex=5,ylim=c(0,(platewide_ave_vic_pos_FI+(6*sdplatewide_ave_vic_pos_FI))),ylab="000_average vic FI, pos/neg and 99.9% CI",xaxt="n")
points(samples$ave_fi_vic_negatives,pch=46,cex=5)
abline(h=(platewide_ave_vic_pos_FI+2.59*sdplatewide_ave_vic_pos_FI), col="PINK", lty=1)
abline(h=(platewide_ave_vic_pos_FI-2.59*sdplatewide_ave_vic_pos_FI), col="PINK", lty=1)
abline(h=(platewide_ave_vic_neg_FI+2.59*sdplatewide_ave_vic_neg_FI), col="PINK", lty=1)
abline(h=(platewide_ave_vic_neg_FI-2.59*sdplatewide_ave_vic_neg_FI), col="PINK", lty=1)
axis(1,rownumbers,labels=samples$Var1,cex.axis=0.4,las=2)
grid(nx=100,ny=NULL)
dev.off()

#plot fam average FI for positives and negatives with platewide 99%CI

pdf(file.path(outpath,"000_average_fam_FI.pdf"))

plot(samples$ave_fi_fam_positives,pch=46,col="red",cex=5,ylim=c(0,(platewide_ave_fam_pos_FI+1000)),ylab="average fam FI, pos/neg and 99% CI",xaxt="n")
points(samples$ave_fi_fam_negatives,pch=46,cex=5)
abline(h=(platewide_ave_fam_neg_FI+2.59*sdplatewide_ave_fam_neg_FI), col="PINK", lty=1)
abline(h=(platewide_ave_fam_neg_FI-2.59*sdplatewide_ave_fam_neg_FI), col="PINK", lty=1)
abline(h=(platewide_ave_fam_pos_FI+2.59*sdplatewide_ave_fam_pos_FI), col="PINK", lty=1)
abline(h=(platewide_ave_fam_pos_FI-2.59*sdplatewide_ave_fam_pos_FI), col="PINK", lty=1)
axis(1,rownumbers,labels=samples$Var1,cex.axis=0.4,las=2)
grid(nx=100,ny=NULL)
dev.off()


#assign warnings based on FI average values and 95% CI platewide
samples$vicwarningspositives<-NA
for (i in rownumbers) {if (is.na(samples$ave_fi_vic_positives[i])){next}else{if (samples$ave_fi_vic_positives[i]>platewide_ave_vic_pos_FI+(2.59*sdplatewide_ave_vic_pos_FI)){samples$vicwarningspositives[i]<-"vic average FI positive value is greater than 99% upper CI limit for the plate"}}}
for (i in rownumbers) {if (is.na(samples$ave_fi_vic_positives[i])){next}else{if(samples$ave_fi_vic_positives[i]<platewide_ave_vic_pos_FI-(2.59*sdplatewide_ave_vic_pos_FI)){samples$vicwarningspositives[i]<-"vic average FI positive value is lesser than 99% upper CI limit for the plate"}}}
for (i in rownumbers) {if (is.na(samples$ave_fi_vic_positives[i])){samples$vicwarningspositives[i]<-"ENDOGENOUS CONTROL FAILED"}}



samples$vicwarningsnegatives<-NA
for (i in rownumbers) {if (is.na(samples$ave_fi_vic_negatives[i])){next}else{if (samples$ave_fi_vic_negatives[i]>platewide_ave_vic_neg_FI+(2.59*sdplatewide_ave_vic_neg_FI)){samples$vicwarningsnegatives[i]<-"vic average FI negative value is greater than 99% upper CI limit for the plate"}}}
for (i in rownumbers) {if (is.na(samples$ave_fi_vic_negatives[i])){next}else{if (samples$ave_fi_vic_negatives[i]<platewide_ave_vic_neg_FI-(2.59*sdplatewide_ave_vic_neg_FI)){samples$vicwarningsnegatives[i]<-"vic average FI negative value is lesser than 99% upper CI limit for the plate"}}}
for (i in rownumbers) {if (is.na(samples$ave_fi_vic_positives[i])){samples$vicwarningsnegatives[i]<-"NO BASELINE : DROPLET SATURATION, SYSTEM ERROR OR NO SAMPLE"}}


samples$famwarningspositives<-NA
for (i in rownumbers) {if (is.na(samples$ave_fi_fam_positives[i])){next}else{if(samples$ave_fi_fam_positives[i]>platewide_ave_fam_pos_FI+(2.59*sdplatewide_ave_fam_pos_FI)){samples$famwarningspositives[i]<-"fam average FI positive value is greater than 99% upper CI limit for the plate"}}}
for (i in rownumbers) {if (is.na(samples$ave_fi_fam_positives[i])){next}else{if(samples$ave_fi_fam_positives[i]<platewide_ave_fam_pos_FI-(2.59*sdplatewide_ave_fam_pos_FI)){samples$famwarningspositives[i]<-"fam average FI positive value is lesser than 99% upper CI limit for the plate"}}}

samples$famwarningsnegatives<-NA
for (i in rownumbers) {if (is.na(samples$ave_fi_fam_negatives[i])){next}else{if (samples$ave_fi_fam_negatives[i]>platewide_ave_fam_neg_FI+(2.59*sdplatewide_ave_fam_neg_FI)){samples$famwarningsnegatives[i]<-"fam average FI negative value is greater than 99% upper CI limit for the plate"}}}
for (i in rownumbers) {if (is.na(samples$ave_fi_fam_negatives[i])){next}else{if (samples$ave_fi_fam_negatives[i]<platewide_ave_fam_neg_FI-(2.59*sdplatewide_ave_fam_neg_FI)){samples$famwarningsnegatives[i]<-"fam average FI negative value is lesser than 99% upper CI limit for the plate"}}}
for (i in rownumbers) {if (is.na(samples$ave_fi_fam_negatives[i])){samples$famwarningsnegatives[i]<-"NO BASELINE : DROPLET SATURATION, SYSTEM ERROR OR NO SAMPLE"}}


#PLOT BEAD COUNTS FOR ALL SAMPLES
platewide_average_beadcount<-mean(samples$Freq)
platewide_sd_beadcount<-sd(samples$Freq)
pdf(file.path(outpath,"000_average_beadcounts.pdf"))
plot(samples$Freq,pch=46,col="red",cex=5,ylim=c(0,(platewide_average_beadcount+(6*platewide_sd_beadcount))),ylab="average droplet count and 99% CI",xaxt="n")
abline(h=(platewide_average_beadcount+2.59*platewide_sd_beadcount), col="PINK", lty=1)
abline(h=(platewide_average_beadcount-2.59*platewide_sd_beadcount), col="PINK", lty=1)
axis(1,rownumbers,labels=samples$Var1,cex.axis=0.4,las=2)
grid(nx=100,ny=NULL)
dev.off()
#QC warning about bead count
samples$dropletcountwarning<-NA
for (i in rownumbers) {if (samples$Freq[i]>platewide_average_beadcount+(2.59*platewide_sd_beadcount)){samples$dropletcountwarning[i]<-"Droplet count value is greater than 99% upper CI limit for the plate"}}
for (i in rownumbers) {if (samples$Freq[i]<platewide_average_beadcount-(2.59*platewide_sd_beadcount)){samples$dropletcountwarning[i]<-"Droplet count value is lesser than 99% upper CI limit for the plate"}}
for (i in rownumbers) {if (samples$Freq[i]>platewide_average_beadcount+(2.59*platewide_sd_beadcount)){samples$FAILS[i]<-"Droplet count value is greater than 99% upper CI limit for the plate"}}
for (i in rownumbers) {if (samples$Freq[i]<platewide_average_beadcount-(2.59*platewide_sd_beadcount)){samples$FAILS[i]<-"Droplet count value is lesser than 99% upper CI limit for the plate"}}
for (i in rownumbers) {if (samples$Freq[i]<10000){samples$FAILS[i]<-"Droplet count is below 10000"}}

rm(current_na)

current_na<-is.na(samples$Freq)
for (i in rownumbers) {if (current_na[i]==TRUE){samples$dropletcountwarning[i]<-"No droplets"}}
for (i in rownumbers) {if (current_na[i]==TRUE){samples$FAILS[i]<-"No droplets"}}


#calculate values for poisson calculation.
#phat is estimator of P, the probability of a micelle being positive. P is unknown but phat can be estimated by phat = number of positive droplets / total droplets
samples$phat_fam<-((samples$famposvicneg+samples$famposvicpos)/samples$Freq)
samples$phat_vic<-((samples$famnegvicpos+samples$famposvicpos)/samples$Freq)

#standard deviation of p_hat estimator
samples$sd_phat_fam<- sqrt((samples$phat_fam*(1-samples$phat_fam))/samples$Freq)
samples$sd_phat_vic<- sqrt((samples$phat_vic*(1-samples$phat_vic))/samples$Freq)

#upper and lower confidence intervals (95%) of p_hat estimates
samples$phat_low_fam <- samples$phat_fam - (1.96 * samples$sd_phat_fam)
samples$phat_high_fam <- samples$phat_fam + (1.96 * samples$sd_phat_fam)

samples$phat_low_vic <- samples$phat_vic - (1.96 * samples$sd_phat_vic)
samples$phat_high_vic <- samples$phat_vic + (1.96 * samples$sd_phat_vic)

#lambda is the true concentration of target molecules per chamber. It is unknown but can be estimated by lambda hat. lambda_hat =  -ln (1-p_hat)

samples$lambda_hat_fam<-(-log(1-samples$phat_fam))
samples$lambda_hat_vic<-(-log(1-samples$phat_vic))

#upper and lower confidence intervals for lambda hat :
samples$lambda_hat_low_fam <- (-log(1-samples$phat_low_fam))
samples$lambda_hat_high_fam <- (-log(1-samples$phat_high_fam))

samples$lambda_hat_low_vic <- (-log(1-samples$phat_low_vic))
samples$lambda_hat_high_vic <- (-log(1-samples$phat_high_vic))

# lambda_hat = concentration target per micelle
# lambda_hat_low = lower 95% CI for lambda hat
# lambda_hat_high = upper 95% CI for lambda hat

# calculate concentration of target molecules per reaction.

samples$fam_copies_per_reaction<-samples$lambda_hat_fam*samples$Freq
samples$vic_copies_per_reaction<-samples$lambda_hat_vic*samples$Freq

samples$fam_copies_per_reaction_low<-samples$lambda_hat_low_fam*samples$Freq
samples$fam_copies_per_reaction_high<-samples$lambda_hat_high_fam*samples$Freq


samples$vic_copies_per_reaction_low<-samples$lambda_hat_low_vic*samples$Freq
samples$vic_copies_per_reaction_high<-samples$lambda_hat_high_vic*samples$Freq




samples$volume<-samples$Freq*0.91e-3
samples$concentration_fam_average_copies_uL <- samples$fam_copies_per_reaction/samples$volume
samples$concentration_fam_low_copies_uL <- samples$fam_copies_per_reaction_low/samples$volume
samples$concentration_fam_high_copies_uL <- samples$fam_copies_per_reaction_high/samples$volume

samples$concentration_vic_average_copies_uL <- samples$vic_copies_per_reaction/samples$volume
samples$concentration_vic_low_copies_uL <- samples$vic_copies_per_reaction_low/samples$volume
samples$concentration_vic_high_copies_uL <- samples$vic_copies_per_reaction_high/samples$volume

for (i in rownumbers){if(is.na(samples$concentration_fam_high_copies_uL[i])){samples$concentration_fam_high_copies_uL[i]<-samples$concentration_fam_average_copies_uL[i]}else{next}}
for (i in rownumbers){if(is.na(samples$concentration_vic_high_copies_uL[i])){samples$concentration_vic_high_copies_uL[i]<-samples$concentration_vic_average_copies_uL[i]}else{next}}

for (i in rownumbers){if(is.infinite(samples$concentration_fam_average_copies_uL[i])){samples$concentration_fam_average_copies_uL[i]<-20000}}
for (i in rownumbers){if(is.infinite(samples$concentration_fam_low_copies_uL[i])){samples$concentration_fam_low_copies_uL[i]<-20000}}
for (i in rownumbers){if(is.infinite(samples$concentration_fam_high_copies_uL[i])){samples$concentration_fam_high_copies_uL[i]<-20000}}
for(i in rownumbers){if(is.infinite(samples$concentration_vic_average_copies_uL[i])){samples$concentration_vic_average_copies_uL[i]<-20000}}
for(i in rownumbers){if(is.infinite(samples$concentration_vic_low_copies_uL[i])){samples$concentration_vic_low_copies_uL[i]<-20000}}
for(i in rownumbers){if(is.infinite(samples$concentration_vic_high_copies_uL[i])){samples$concentration_vic_high_copies_uL[i]<-20000}}



#plot graphs of calculated concenentration (copies/uL) values with error bars
#library(gplots)

pdf(file.path(outpath,"001_concentration_vic.pdf"))
plot (samples$concentration_vic_average_copies_uL,pch=46,cex=4,log="y",xlab="Sample",ylab="HURNASE (copies/uL)",xaxt="n",cex.axis=0.4,xlim=c(0,100),ylim=c(0.1,(max(samples$concentration_vic_high_copies_uL)*1.1)))
for(i in rownumbers){arrows(i,samples$concentration_vic_average_copies_uL[i],i,samples$concentration_vic_high_copies_uL[i],length=0.02,angle=90,code=2,col="red")}
for(i in rownumbers){arrows(i,samples$concentration_vic_average_copies_uL[i],i,samples$concentration_vic_low_copies_uL[i],length=0.02,angle=90,code=2,col="red")}
axis(1,rownumbers,labels=samples$Var1,cex.axis=0.4)
grid()
dev.off()




pdf(file.path(outpath,"001_concentration_fam.pdf"))
plot (samples$concentration_fam_average_copies_uL,pch=46,cex=4,log="y",xlab="Sample",ylab="HURNASE (copies/uL)",xaxt="n",cex.axis=0.4,xlim=c(0,100),ylim=c(0.1,(max(samples$concentration_fam_high_copies_uL)*1.1)))
for(i in rownumbers){arrows(i,samples$concentration_fam_average_copies_uL[i],i,samples$concentration_fam_high_copies_uL[i],length=0.02,angle=90,code=2,col="red")}
for(i in rownumbers){arrows(i,samples$concentration_fam_average_copies_uL[i],i,samples$concentration_fam_low_copies_uL[i],length=0.02,angle=90,code=2,col="red")}
axis(1,rownumbers,labels=samples$Var1,cex.axis=0.4)
grid()
dev.off()



#DEFINE FAILS AND RESULTS
for (i in rownumbers) {if (samples$concentration_fam_low_copies_uL[i]<0){samples$FAILS[i]<-"WARNING : CI for FAM crosses zero"}}
for (i in rownumbers) {if (samples$famnegvicneg[i]<(0.01*samples$Freq[i])){samples$FAILS[i]<-"WARNING: REACTION MAY BE SATURATED WITH TOO MANY TEMPLATES"}}
for (i in rownumbers) {if (samples$concentration_vic_low_copies_uL[i]<0.0000001){samples$FAILS[i]<-"FAIL: HUMAN DNA CONCENTRATION TOO LOW"}}

#define amount of area under the curve that is above zero
for (i in rownumbers){samples$areaabovezerofam[i]<-1-pnorm(0,mean=samples$phat_fam[i],sd=samples$sd_phat_fam[i])}
for (i in rownumbers){samples$areaabovezerovic[i]<-1-pnorm(0,mean=samples$phat_vic[i],sd=samples$sd_phat_vic[i])}

#assign endogenous control status
samples$endogenouscontrolresult<-NA
for (i in rownumbers) {if (samples$areaabovezerovic[i]>=0.95){samples$endogenouscontrolresult[i]<-"Endogenous control OK"}}
for (i in rownumbers) {if (samples$areaabovezerovic[i]<0.90){samples$endogenouscontrolresult[i]<-"Endogenous control failed"}}

#assign warnings based on FI average values and 95% CI platewide
samples$targetresult<-NA
for (i in rownumbers) {if (samples$areaabovezerofam[i]>=0.95){samples$targetresult[i]<-"Sample is positive for the target"}}
for (i in rownumbers) {if (samples$areaabovezerofam[i]<0.90){samples$targetresult[i]<-"Sample is negative for the target"}}




#print output

write.table(samples, file = (file.path(outpath,"000_results_full.txt")),row.names=F,col.names=T,quote=F,sep="\t")
columns<-c(1,55,56,47,48,49,50,51,52,53,54,2,13,14,15,16,11,12,46)
write.table(samples[,columns], file = (file.path(outpath,"000_results_concise.txt")),row.names=F,col.names=T,quote=F,sep="\t")


##########################################################################################
#EXPLANATION OF CONCISE REPORT
#
#Var1 									: The name of the specimen (plate number and well location)
#endogenouscontrolresult 				: Was the endogenous control positive (VIC zeta score at least 0.95)?
#targetresult 							: Was the sample positive for the target (FAM zeta at least 0.95)?
#concentration_fam_average_copies_uL	: Mean average estimated concentration target (copies/uL of the PCR mix)
#concentration_fam_low_copies_uL		: Lower 95% CI estimated concentration target  (copies/uL of the PCR mix)
#concentration_fam_high_copies_uL		: Upper 95% CI estimated concentration target  (copies/uL of the PCR mix)
#concentration_vic_average_copies_uL	: Mean average estimated concentration endogenous control (copies/uL of the PCR mix)
#concentration_vic_low_copies_uL		: Lower 95% CI estimated concentration endogenous control (copies/uL of the PCR mix)
#concentration_vic_high_copies_uL		: Upper 95% CI estimated concentration endogenous control (copies/uL of the PCR mix)
#areaabovezerofam						: Zeta value (FAM)
#areaabovezerovic						: Zeta value (VIC)
#Freq									: Total number of assayed droplets
#volume									: Total volume of PCR mixture assayed
#famposvicneg							: Count of FAM positive, VIC negative droplets
#famnegvicpos							: Count of FAM negative, VIC positive droplets
#famposvicpos							: Count of FAM positive, VIC positive droplets	
#famnegvicneg							: Count of FAM negative, VIC negative droplets	
#famthreshold							: Threshold value of FAM fluorescence intensity for positive droplet classification
#victhreshold							: Threshold value of VIC fluorescence intensity for positive droplet classification	
#
#The full report contains more detailed information and shows how copies/uL eastimates relate to the original specimen.
#In most cases, it will not be necessary to refer to the full report.
##########################################################################################
#END OF SCRIPT
##########################################################################################

