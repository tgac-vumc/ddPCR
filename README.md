# README file for ddPCR
# 
# 201605, HFvanEssen

Use ddPCR QuantaSoft software to export raw Amplitude data
- 'Export Amplitude and Cluster Data' in options 

Pipeline example:
- createDesign
- loop through individual targets
- combine all samples for target 
- find outliers
- find global thresholds combined data
- plot combined data with basic thresholds
- refine global thresholds for combined data
- define clusters for combined data
- plot combined data with refined thresholds
- loop through each sample for each target with refined thresholds
- plot data
- use calculateStats... functions to get sample statistics
######