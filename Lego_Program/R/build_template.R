build_template <- function(num_players, pw_prob, tr.avoid) {
  
  
  
  #Define player interaction matrix
  int_m <- matrix(0,num_play,num_play)
  
  
  #Fill out matrix based on probabilities above
  prob_line <- cumsum(sort(pw_prob))
  
  for (i in 1:num_play) {
    for (j in 1:num_play) {
      
      if (i < j) {
        
        #draw an interaction pair
        r_draw <- runif(1)
        
        if (r_draw < prob_line[2]) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[1],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[1]) && (r_draw < prob_line[2])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[2],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[2]) && (r_draw < prob_line[3])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[3],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[3]) && (r_draw < prob_line[4])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[4],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[4]) && (r_draw < prob_line[5])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[5],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[5]) && (r_draw < prob_line[6])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[6],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[6]) && (r_draw < prob_line[7])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[7],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[7]) && (r_draw < prob_line[8])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[8],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
        if ((r_draw > prob_line[8]) && (r_draw < prob_line[9])) {
          #Obtain the interaction types from the probability line, randomize, and assign
          int_name <- unlist(strsplit(strsplit(names(prob_line)[9],"_")[[1]][2],""))
          mij <- sample(int_name,2,replace=FALSE)
          int_m[i,j] <- mij[1]
          int_m[j,i] <- mij[2]
        }
        
      }
    }
  }
  #Determine which are active players and which aren't
  for (i in 1:num_play) {
    if (length(which(int_m[i,] == "e")) > 0) {
      int_m[i,i] <- "n"
    } else {
      int_m[i,i] = "i"
    }
  }
  
  #Implement posthoc Rules
  
  #1: Row/Col 1 is the sun
  int_m[1,] <- rep("i",num_play)
  #Things can only eat the sun... they cannot avoid, need, or make it
  #suggestion:
  int_m[which(int_m[,1]!="e"),1]="i"
  
  #2: if player A contains an "m" with player B, player B is "i" with everything
  # except it is "n" with player A
  #Which species 'make things?'
  #NOTE: multiple A,B,C can make the same D
  
  for (i in 1:num_play) {
    #what things are made?
    made <- which(int_m[i,] == "m")
    l_made <- length(made)
    if (l_made > 0) {
      for (k in 1:l_made) {
        int_m[made[k],] <- rep("i",num_play) #ignores all
        #what thing(s) make it?
        makers <- which(int_m[,made[k]] == "m")
        int_m[made[k],makers] <- "n" #Except the thing that makes it
      }
    }
  }
  
  #Determine ACTIVE player similarity for avoidance
  
  labels <- unlist(sapply(seq(1,num_play),function(x){paste("P",x,sep="")}))
  colnames(int_m) <- labels
  rownames(int_m) <- labels
  
  #One each time so avoidance can be asymmetrical; 
  #tr.avoid=0.8 #the threshold determines how much overlap is needed for avoidance
  #This seems to be very similar to tilman rules of competition:
  #... a species excludes the other if it exploits (partially) the resources that the other needs the most
  for (k in 2:num_play){ #Start from column 2 to avoid competition for the sun
    for (l in 3:(num_play-1)){
      aux=int_m[k,] #fix a row
      aux.e=which(aux=="e") #positions of resources
      aux.n=which(aux=="n") #position of needs
      
      aux.c=int_m[l,c(aux.e,aux.n)] #the potential "competitor"
      sum.rn=length(c(aux.e,aux.n)) #number of resources and needs
      
      if (sum.rn>0){ #avoiding div by 0 (things may not have needs or eat)
        share=sum(aux[c(aux.e,aux.n)]==aux.c) #number of shared resources and needs
        prop.sh=share/sum.rn
        if (prop.sh > tr.avoid){ #only avoids if proportion of shared needs/rec > threshold
          int_m[k,l]="a"
        }  
      }
    }
  }
  
  #Determine which are active players and which aren't
  for (i in 1:num_play) {
    if (length(which(int_m[i,] == "e")) > 0) {
      int_m[i,i] <- "n"
    } else {
      int_m[i,i] = "i"
    }
  }
  
  
  return(int_m)
  
  
  
}