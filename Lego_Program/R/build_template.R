build_template <- function(num_players, pw_prob) {
  
  
  
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
  
  return(int_m)
  
  
  
}