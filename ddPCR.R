# PACKAGES

if("dplyr" %in% installed.packages()[,1] == FALSE){
  install.packages("dplyer")
} else { library(dplyr) }
if("magrittr" %in% installed.packages()[,1] == FALSE){
  install.packages("magrittr")
} else { library(magrittr) }
if("xlsx" %in% installed.packages()[,1] == FALSE){
  install.packages("xlsx")
} else { library(xlsx) }

#library(ddpcR)
# end, HF van Essen 2016
  