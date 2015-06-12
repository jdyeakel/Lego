test_compatability <- function(orig_subint_m,labels) {
  
  subint_m <- orig_subint_m
  compatable <- FALSE
  
  while (compatable == FALSE) 
    {
    
    num_play <- length(subint_m[,1])
    
    for (i in 2:num_play)  #ignore the sun
      {
      #Test compatability of PASSIVE players
      if (subint_m[i,i] == "i") 
        {
        #Where are things making the passive player?
        makers <- labels[which(subint[i,] == "n")]
        #Are the makers present?
        test <- any(makers %in% colnames(subint_m))
        if (test == FALSE) 
          {
          subint_m <- subint_m[-i,]
          subint_m <- subint_m[,-i]
          break
        }
      } else #Test compatability of ACTIVE players
        {
          
          #Does the species eat something?
          if (length(which(subint_m[i,] == "e")) == 0) {
            subint_m <- subint_m[-i,]
            subint_m <- subint_m[,-i]
            break
          } 
          #Does the species have what it needs?
          required <- labels[which(subint_m[i,] == "n")]
          if (length(required) > 0) {
            test <- all(required %in% colnames(subint_m))
            if (test == FALSE) 
            {
              subint_m <- subint_m[-i,]
              subint_m <- subint_m[,-i]
              break
            }
          }
      }
      
    } #End For Loop
    
    if (i == num_play) {
      compatable <- TRUE
    } else {
      compatable <- FALSE
    }
    
  } #End While loop
  
  return subint_m
  
} #End function