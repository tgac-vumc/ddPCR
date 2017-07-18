### density lines ------------------------------------------

# read in data of of one sample 
# determine the 


kde <- kde2d(temp$nReads[finished],temp$ClusterDensity[finished], n = 75,lims=c(-20,40,-200,1700))
contour(kde$x, kde$y, kde$z, add = TRUE, col = "red", lwd = 2, drawlabels = FALSE)
points(temp$nReads[!finished],temp$ClusterDensity[!finished],pch=20,col="green",cex=1.4)

### slope or angle  ------------------------------------------

