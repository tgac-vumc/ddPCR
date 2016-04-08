createTargetFolders <- function(x, path=NULL)
{ # x = output of get.targets(). Will create output folders at specified path.
  paths <- NULL
  if(class(x) == "list" & class(path) != "NULL")
  {
    paths <- file.path(path,names(x))
    for(i in 1:length(paths))
    {
      if(!file.exists(paths[i]))
      {
        dir.create(paths[i])
      }
    }
  } else if (class(x) == "character" | class(x) == "factor")
  {
    if(class(path) != "NULL")
    {
      for(i in 1:length(x))
      {
        paths[i] <- file.path(path,x[i])
        if(!file.exists(paths[i]))
        {
          dir.create(paths[i])
        }
      }
    }
  }
  return(paths)
}