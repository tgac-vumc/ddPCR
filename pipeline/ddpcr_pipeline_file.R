# This is and example of the ddPCR pipeline analysis.
# Before use, the file MUST be copied to a different
# location. Then the "#" may be removed and the code
# can be run. Please do not make any changes in this
# document.
#
# HFvE, 201711

#### USER INPUT OSX
#source.path <- "/Users/dirkvanessen/Documents/coding/R/ddPCR-s4/ddPCR"
#data.path <- "/Users/dirkvanessen/Desktop/20170725 EGFR ddPCR L858R en T790M MvM en MD_2017-07-25-15-56"

#### USER INPUT WINDOWS
source.path <-"\\\\vumc.nl/onderzoek$/s4e-gpfs1/pa-tgac-01/analisten/Dirk/r_scripts/ddpcr_analysis/scripts_ddpcr"
data.path <- "\\\\vumc.nl/home$/store4ever/h.vanessen/h.vanessen/Desktop/20171106 EGFRs ddPCRs SV_2017-11-06-16-18"

#### SOURCE FUNCTION
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    if(trace) cat(nm,"")           
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir(source.path)

### READING AMPLITUDE FILE
files <- list.files(data.path, pattern = "Amplitude.csv", full.names = TRUE)

data <- readAmplitudeFiles(files, verbose = TRUE)
data <- setThresholds(data = data, algorithm = "densityhist", 
                        breaks = 25, strict = TRUE, rm.outliers = TRUE,
                        verbose = TRUE)

overview(data)

#data <- minOutliers(data, breaks = 300, strict = TRUE)

plot.ddPCR(data, well = 1, new = TRUE, verbose = TRUE)
plot.ddPCR(data, well = 3, new = TRUE, verbose = TRUE)
plot.ddPCR(data, well = 5, new = TRUE, verbose = TRUE)

overview(data)

### show all plots
for(i in 1:12{
  plot.ddPCR(data, well = i, new = TRUE)
  Sys.sleep(0.5)
}

for(i in 1:length(files)){
  output.file <- gsub(basename(files[i]), pattern = "_Amplitude.csv", replacement = ".png")
  outputFile <- file.path(data.path, output.file)
  png(filename = outputFile, width = 750, height = 750)
  plot.ddPCR(data, well = i, new = TRUE)
  dev.off()
}
### 
saveRDS(data, file = "ddPCR_data_file.Rds")
test <- readRDS(file = "ddPCR_data_file.Rds")