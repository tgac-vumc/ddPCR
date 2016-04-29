plotThresholds <- function(thresholds, col = "black")
{
  abline(h = thresholds[1], col = col)
  abline(v = thresholds[2], col = col)
}