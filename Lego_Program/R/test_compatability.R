test_compatability <- function(orig_subint_m,r_sample) {
  
  subint_m <- orig_subint_m[r_sample,r_sample]
  labels<-row.names(subint_m)
  
  compatable <- FALSE

  
  while (compatable == FALSE) 
    {
    
    num_play <- length(subint_m[,1])
    
    for (i in 2:num_play)#ignore the sun
      {  
      orig.pos=which(row.names(orig_subint_m)==labels[i]) #position in the original matrix
      
      #Test compatability of PASSIVE players
      if (subint_m[i,i] == "i") 
        {
        #is it made by someone?
        if (sum(subint_m[i,] == "n")>0){ #previously it was considering the below logical op even if there were no makers
          #Where are things making the passive player (we have always to refer to the original)?
          makers <- names(which(orig_subint_m[orig.pos,] == "n"))
          #Are the makers present?
          test <- any(makers %in% colnames(subint_m)) 
          if (test == FALSE) 
          {
            if (length(subint_m)==4) {break; compatable=TRUE} #avoids error with 2x2 matrix
            subint_m <- subint_m[-i,]
            subint_m <- subint_m[,-i]
            break
          }
        }
      } else #Test compatability of ACTIVE players
      {
        
        #Does the species eat something?
        if (length(which(subint_m[i,] == "e")) == 0) {
          if (length(subint_m)==4) {break; compatable=TRUE} #avoids error with 2x2 matrix
          subint_m <- subint_m[-i,]
          subint_m <- subint_m[,-i]
          break
        } 
        #Does the species have what it needs (we have always to refer to the original)?
        required <- names(which(orig_subint_m[orig.pos,] == "n"))
        if (length(required) > 1) { #>1 because it always need itself
          test <- all(required %in% colnames(subint_m))
          if (test == FALSE) 
          {
            if (length(subint_m)==4) {break; compatable=TRUE} #avoids error with 2x2 matrix
            subint_m <- subint_m[-i,]
            subint_m <- subint_m[,-i]
            break
          }
        }
      }
    } #End For Loop
    
    if (i == num_play) {
      compatable <- TRUE
    #} else {
     # compatable <- FALSE
    }
    
  } #End While loop
  
  return (subint_m)
  
} #End function