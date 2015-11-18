# README file for ddPCR
# https://github.com/drkvnssn/ddPCR-pipeline.git
# 201511, HFvanEssen

Files needed to run an update on MiSeq output folder

- Use ddPCR QuantaSoft software to analyse the samples
- Export as csv
-- What is native? Is it usable?
- csv per sample well
	- contains amplitude plus cluster definitions.
- Create design.file for experiment
	 [user created, contains: sample name, filename, sample type, and probe].
- Read in csv together with design.file
- One design file to catalogue all csv files.
- Use ddPCR scripts to analyse the data.
- Pos and neg controls are combined and analysed. - determines breakpoints for clusters
- Breakpoints are used to analyse the samples.
- Plots!


