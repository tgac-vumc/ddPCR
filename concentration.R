#function concentration(pos_count,total_count){
#  var ncount = total_count - pos_count;
#  +	// this is the key calculation whic returns the copy number per 20 micro litres.
#  var ratio = Math.abs(-Math.log(ncount/total_count)/0.91)*1000*20; 
#  //console.debug(ratio);
#  return ratio; 
  
  
concentration <-  function(negCount, Count, vDroplet=0.91, volume=1)
  { # concentration in copies / user defined volume
    result <- ((-log(negCount/Count)/vDroplet))*1000*volume
    return(result)
  }