createDesign <- function(path, probe = "probe")
{
  amplitude.files <- list.files(path = path, pattern = "_Amplitude.csv")
  sample.names <- gsub(pattern = "_Amplitude.csv", x = amplitude.files, replacement = "")
  design <- matrix(data = "", nrow = length(amplitude.files), ncol = 4)
  colnames(design) <- c("Name", "File", "Type", "Probe")
  design[,1] <- sample.names
  design[,2] <- amplitude.files
  design[,3] <- c("pos", (rep("sample", length(amplitude.files)-2)), "neg")
  design[,4] <- probe
  output.file <- file.path(path, "design.txt")
  write.table(file = output.file, x = design, quote = FALSE, sep = "\t", row.names = FALSE)
}