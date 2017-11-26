overview <- function(data){
  if((class(data)[1] == "ddPCRdata") != TRUE){
    stop ("data structure is not in the correct format.\n")
  } 
  probes <- unique(data@phenoData$sampleData['probe', ])
  nSamples <- ncol(data@phenoData$sampleData)
  cat("ddPCR structure has", nSamples, "samples with", length(probes), "probe(s).\n\n")
  for (i in 1:length(probes)){
    selection <- data@phenoData$sampleData['probe', ] %in% probes[i]
    message("Probe ",i, ": ") ; cat(probes[i], "\n\n", sep ="")
    
    rownames <- c(1:nSamples)[selection]
    samples <- data@phenoData$sampleData['name', selection ]
    wells <- data@phenoData$sampleData['well', selection ]
    threshold.ch1 <- round(data@phenoData$ch1['threshold', selection])
    threshold.ch2 <- round(data@phenoData$ch2['threshold', selection])
    minOutlier.ch1 <- round(data@phenoData$ch1['minOutlier', selection])
    minOutlier.ch2 <- round(data@phenoData$ch2['minOutlier', selection])
    
    results <- matrix(NA, nrow = length(samples), ncol = 6,
                      dimnames = list(rownames,c("Name", "Well", 
                                                  "Thresh.ch1", 
                                                  "Thresh.ch2",
                                                 "minOut.ch1",
                                                 "minOut.ch2")))
    results[,'Name'] <- samples
    results[,'Well'] <- wells
    results[,'Thresh.ch1'] <- threshold.ch1
    results[,'Thresh.ch2'] <- threshold.ch2
    results[,'minOut.ch1'] <- minOutlier.ch1
    results[,'minOut.ch2'] <- minOutlier.ch2
    results <- data.frame(results)
    print(results)
    cat("\n")
  }
}
