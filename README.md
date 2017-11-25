# README file for ddPCR master-S4
# 
# 201711, HFvanEssen
#
# Current branch is a work in progress. Amplitude files from one experiment (plate) can be read. All data is stored into one 'ddPCR' S4 structure. From this structure multiple analysis can be perfomed. 

Use ddPCR QuantaSoft software to export raw Amplitude data
- 'Export Amplitude and Cluster Data' in options 

### select a path and get all the files 
path <- "experiment path"
files <- list.files(path, pattern = "Amplitude.csv", full.names = TRUE)

### Pipeline example:

data <- readAmplitudeFiles(files)

### set a threshold for the samples. Current algorithms: 'histogram', 'kmean2', 'ranges', 'densityhist'
### for algorithms with histogram 
data <- setThresholds(data = data, algorithm = "densityhist", 
                        breaks = break.sizes[z], strict = TRUE, 
                        verbose = FALSE)
data <- minOutliers(data)
plot.ddPCR(data, well = i, new = TRUE)

}
- createDesign
- loop through individual targets
- combine all samples in each used probe
- find outliers
- find global thresholds combined data
- plot combined data with basic thresholds
- refine global thresholds for combined data
- plot combined data with refined thresholds
- loop through each sample for each target with refined threshold
- plot data
- use calculateStats... functions to get sample statistics
######
