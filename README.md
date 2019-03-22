### README file for ddPCR master-S4
201711, HFvanEssen

Current branch is a work in progress. Amplitude files from one experiment (plate) can be read. All data is stored into one 'ddPCR' S4 structure. From this structure multiple analysis can be perfomed. 

# EXPORT AMPLITUDE FILES
Use ddPCR QuantaSoft software to export raw Amplitude data
- 'Export Amplitude and Cluster Data' in options 

# SELECT ALL THE FILES

select a path and get all the files 

```
path <- "experiment path"
files <- list.files(path, pattern = "Amplitude.csv", full.names = TRUE)
```

# READING DATA

```
data <- readAmplitudeFiles(files)
```

If the experiment file (also .csv) is present in the same folder then additional metadata will be imported as well into the S4 structure. Data such as probe names will be used to group similar samples together during analysis.

# Set a threshold for the data
Set a threshold for the samples. Current algorithms: 'histogram', 'kmean2', 'ranges', 'densityhist'.

For algorithms that use histogram as a base the number of breaks can be added and if breaks are strict (that is exactly as predefined by the user). 

Strict = TRUE or FALSE

Example to set thresholds.

``` 
data <- setThresholds(data = data, algorithm = "densityhist", 
                        breaks = break.sizes[z], strict = TRUE, 
                        verbose = FALSE)
```
                
# PLOTTING

After setting the threshold for each sample. The data can be plotted with the function below. Two variations of the plot functions are available (new = TRUE or FALSE)

```
plot.ddPCR(data, well = i, new = TRUE)
```

# PLOT COLOURS

Standard ddpcr colors are used for plotting, but can be changed with changeColors and a string of 6 hex colors in this order

cluster 0 : outlier
cluster 1 : double negative
cluster 2 : ch 1 negative, ch 2 positive
cluster 3 : double positive
cluster 4 : ch 1 postive, ch 2 negative
cluster 5 : rain
