# README file for ddPCR master-S4
# 
# 201711, HFvanEssen
#
# Current branch is a work in progress. Amplitude files from one experiment (plate) can be read. All data is stored into one 'ddPCR' S4 structure. From this structure multiple analysis can be perfomed. 

Use ddPCR QuantaSoft software to export raw Amplitude data
- 'Export Amplitude and Cluster Data' in options 

### select a path and get all the files 
'''
path <- "experiment path"
files <- list.files(path, pattern = "Amplitude.csv", full.names = TRUE)
'''
### Pipeline example:

data <- readAmplitudeFiles(files)

### set a threshold for the samples. Current algorithms: 'histogram', 'kmean2', 'ranges', 'densityhist'
### for algorithms with histogram the number of breaks can be added and if they are to made scrict 
### or not strict (strict = TRUE/FALSE)

data <- setThresholds(data = data, algorithm = "densityhist", 
                        breaks = break.sizes[z], strict = TRUE, 
                        verbose = FALSE)
                        
### After setting the threshold for each sample. The data can be plotted with the function below.
### A new plot function is under development (new = TRUE)

plot.ddPCR(data, well = i, new = TRUE)

### Standard ddpcr colors are used for plotting, but can be changed with changeColors and a string 
### of 6 hex colors in this order
  # cluster 0 : outlier
  # cluster 1 : double negative
  # cluster 2 : ch 1 negative, ch 2 positive
  # cluster 3 : double positive
  # cluster 4 : ch 1 postive, ch 2 negative
  # cluster 5 : rain
  
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
