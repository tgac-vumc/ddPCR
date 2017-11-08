checkPackages <- function(x){
  for (i in 1:length(x)){
    package <- as.character(x[i])
    if( package %in% installed.packages()[,1] == FALSE ){
      install.packages(package, lib = "N:/Documenten/R/win-library/3.2")
    } else { library(package, character.only = TRUE ) }
  }
}