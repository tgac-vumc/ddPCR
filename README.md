# README file for ddPCR
# https://github.com/drkvnssn/ddPCR-pipeline.git
# 201512, HFvanEssen

Files needed to run an update on MiSeq output folder

Use ddPCR QuantaSoft software to analyse the samples
- 'Export Amplitude and Cluster Data' in options 

To run the pipeline select a folder containing:
- experiment setup .csv file
- 'Amplitude.csv'  file for each sample well in the experiment.

Return:
In the experiment path new folders - one for each target - will be created containing:
- ddPCR plot for positive and ntc control.
- ddPCR plot for each sample well in the experiment.
- experiment / target .txt file containing the results of the experiment.



