.removeBelowHighCluster <- function(ch1 = NULL, ch2 = NULL, 
                                 breaks = 100, strict = TRUE,
                                 verbose = FALSE){
  ch1 <- ch1[!(ch1 %in% NA)]
  ch2 <- ch2[!(ch2 %in% NA)]
  if(strict == TRUE){
    breaks.ch1 <- .makeBreaks(min = min(ch1), max = max(ch1), breaks = breaks)
    breaks.ch2 <- .makeBreaks(min = min(ch2), max = max(ch2), breaks = breaks)
  } else {
    breaks.ch1 <- breaks.ch2 <- breaks
  }
  hist.ch1 <- hist(ch1, breaks = breaks.ch1, plot = FALSE)
  high.dens.ch1 <- hist.ch1$mids[hist.ch1$counts == max(hist.ch1$counts)]
  
  hist.ch2 <- hist(ch2, breaks = breaks.ch2, plot = FALSE)
  high.dens.ch2 <- hist.ch2$mids[hist.ch2$counts == max(hist.ch2$counts)]
  
  x <- cbind(ch1, ch2)
  x <- x[x[,1] > high.dens.ch1,]
  x <- x[x[,2] > high.dens.ch2,]
  
  return(x)
}