
# - [ ] create files that are needed for the pipeline

  # "design file" with sample information
    # name  file/well type  probe FDR_positive  FDR_rain
    # name = sample name
    # name/well = sample well name /sample file nae
    # sample type = 'positive', 'sample', 'negative'
    # probe = probe name to select the 'probe file'
    # FDR_positive = the number of positive droplets in the negative control in the 'positive cluster'
    # FDR_rain = the number of droplets in the negative control in the 'rain'

  # probe file with probe information
    # date of each data set
    # breakpoint channel 1, multiple (learning)
    # breakpoint channel 2, multiple (learning)
    # cluster amplitude channel 1, multiple (learning)
    # cluster amplitude channel 2, multiple (learning)
    # cluster sd channel 1, multiple (learning)
    # cluster sd channel 2, multiple (learning)

# - [ ] ideas
  # read all data from all wells and find the breakpoint for all (both channels)
  # there must be 4 clusters for all the data (positive/negative/samples)
  # if auto clustering is done (without set breakpoint), is the breakpoint then at the same level?
  # is the difference between positive amplitude and negative amplitude the same for the control?

# - [ ] function: create basic design file based on files in a folder
  # strsplit and find if they are similar...
  # combine data from multiple plates option?
# - [ ] function: read design
# - [ ] function: check files needed for analysis : sample files, probe file
# - [ ] function: read sample files
# - [ ] function: check positive / negative controls available
# - [ ] function: compare positive / negative controls to probe file data
# - [ ] function: normalize the data based on negative probes ? mean/median/mode option?
# - [ ] function: calculate calls
# - [ ] function: parse probe controls and store data in probe data file
# - [ ] find the min & max for all droplets 
  # use min / max values for plotting of all the wells in this analysis set

# - [ ] combine positive and negative control 
  # calculate calls (see below)
  # calculate the breakpoint
  # calculate the amplitude of cluster 1 & 2, channel 1
  # calculate the amplitude of cluster 1 & 2, channel 2
  # are there 3 or 4 clusters for the positive/negative controls
  # should all samples be normalized based on the mean/median of the negative cluster data?
  # should a flag class be used based on a 2-D plot? 
    # negative droplet channel 1 & negative channel 2
    # negative droplet channel 1 & rain channel 2
    # negative droplet channel 2 & rain channel 1
    # negative droplet channel 1 & positive droplet channel 2
    # negative droplet channel 2 & positive droplet channel 1
    # positive droplet channel 1 & rain channel 2
    # positive droplet channel 2 & rain channel 1
    # positive droplet channel 1 & positive droplet channel 2
    # rain channel 1 & rain channel 2

# - [ ] for each sample
  # find the clusters in the sample based on the positive and negative cluster data and the breakpoint
  # is the positive cluster amplitude at a similar height as the control?
  # calculate the amount of positive droplets -> 1
  # calculate the amount of negative droplets -> -1
  # calculate the amount of droplets rain -> 0
  # add a column to flag the data with the calls (-1,0,1)

# - [ ] calculate the FDR base on the positive and negative controls
  # add calculations/results to probe file 
  # how many of these results should be stored

# - [ ] FDR ANALYSIS
  # read in all positive and negative controls
  # calculate the FDR for all negative samples

# - [ ] 
# - [ ] 
# - [ ] 

