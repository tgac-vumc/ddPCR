# [x] use highest threshold per probe for plotting.
# [X] update clusters after setting thresholds
# [X] update colors after setting clusters
# [x] add function for minOutliers
# [x] update clusters after setting minOutliers
# [x] update colors after setting minOutliers
# [ ] add function for defineTheRain
# [ ] update clusters after setting defineTheRain
# [ ] update colors after setting defineTheRain
# [x] add probe name to the plot
# [x] update overview with correct sample locations as row names
# [ ] add Cluster 1 mean to 'overview' function
# [ ] fix bug 'minOutliers'
# [ ] determine the best amount of breaks for densityhist -> works good >= 25
# [ ] integrate mean + 3sd
# [ ] integrate 'outliers' package


### idea new to remove outliers?
# [ ] mean first big cluster, 
# [ ] mirror everything left to right side = neg cluster
# [ ] take mean + 2.5sd or 3sd to determine the edges
# [ ] all droplets < edge = minOutlier

# or use 'outliers'
library(outliers)
# select negative ch1
# select negative ch2
# outliers in both -> if true in both... remove
