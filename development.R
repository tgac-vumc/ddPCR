get.wells <- function(prefix=NULL,letters=8, numbers=12,suffix=NULL)
{
  if(class(letters) != "numeric" | class(numbers) != "numeric"){
    stop("input values should be numeric.\n")
  } else {
    result <- rep(NA, (letters*numbers))
    location <- c(1,numbers)
    for(i in 1:letters)
    {
      result[location[1]:location[2]] <- paste(LETTERS[i],c(as.character(1:numbers)), sep="")
      location <- location + numbers
    }
  } 
  if(class(prefix) != "NULL")
  {
    result <- paste(prefix,result,sep="")
  }
  if(class(suffix) != "NULL")
  {
    result <- paste(result,suffix,sep="")
  }
  return(result)
}


      
      