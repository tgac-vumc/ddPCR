overview <- function(data){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  } 
  probes <- unique(data@phenoData$sampleData['probe', ])
  nSamples <- ncol(data@phenoData$sampleData)
  cat("ddPCR structure has", nSamples, "samples with", length(probes), "probe(s).\n\n")
  for (i in 1:length(probes)){
    selection <- data@phenoData$sampleData['probe', ] %in% probes[i]
    cat("Probe ",i ,": " , probes[i], "\n", sep ="")
    samples <- data@phenoData$sampleData['name', selection ]
    wells <- data@phenoData$sampleData['well', selection ]
    threshold.ch1 <- data@phenoData$ch1['threshold', selection ]
    threshold.ch2 <- data@phenoData$ch2['threshold', selection ]
    results <- matrix(NA, nrow = length(samples), ncol = 4,
                      dimnames = list(c(1:length(samples)),c("Name", "Well", 
                                                             "Threshold.ch1", "Threshold.ch2")))
    results[,'Name'] <- samples
    results[,'Well'] <- wells
    results[,'Threshold.ch1'] <- threshold.ch1
    results[,'Threshold.ch2'] <- threshold.ch2
    print(results)
    cat("\n")
  }
}
