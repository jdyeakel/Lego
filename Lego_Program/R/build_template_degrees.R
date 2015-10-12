build_template <- function(num_players, pw_prob, tr.avoid) {
  
  num_play <- num_players   
  
  
  #Fill out matrix based on probabilities above
  prob_line <- cumsum(sort(pw_prob)) 
  
  # Distributing trophic interactions (e's) according to a degree distribution
  #probabilities are for the whole matrix
  #this means we should have N.e<-p.e*(num_play*(num_play-1)) trophic interactions (e)
  #thus the mean degree should be mean.k=N.e/num_play
  #we can sample degrre from a exponential ditribution with mean (1/rate) equal to mean.k; degrees=rexp(num_play,1/mean.k)
  
  
  #Still, the degree distrobution scales with fw size.
  N.e<-p.e*(num_play*(num_play-1))
  mean.k<-N.e/(num_play-1)
  aux=1
  while(aux==1){
    degrees<-rexp(num_play-1,1/mean.k)
    degrees<-round(degrees)
    degrees=c(0,degrees) #degree of the sun is 0
    if(max(degrees)<num_play-1){aux=0}
  }
  
  #determining trophic interactions in the interaction matrix
  int_m <- matrix(0,num_play,num_play)
  for (i in 2:num_play) {
    #for (j in 1:num_play) {
    vec<-(1:num_play)
    vec<-vec[-i]
    resource<-sample(vec,degrees[i])
    int_m[i,resource]<-"e"
    
    #Assigning complimentary interactions
    ee_int=which(int_m[resource,i]=="e") # first disconsider ee interactions
    if(length(ee_int)>=1){resource=resource[-ee_int]}
    
    bin.draw<-rbinom(length(resource),1,pw_prob[1]/pw_prob[6]) #randomly determine ei or en (ee already emerge from assigning e's)
    #Faculatative mutualism
    int_m[resource*bin.draw,i]<-"n" #assigning n's according to random draw
    #Asymmetric predation
    int_m[resource*(1-bin.draw),i]<-"i"#assigning i's according to random draw
  }
  
  
  #We now focus on non-trophic interactions
  #First, remove probabilities including 'e'...
  pw_prob.new=pw_prob[-c(1,6,9)]
  pw_prob.new=pw_prob.new/sum(pw_prob.new)
  prob_line <- cumsum(sort(pw_prob.new)) 
  
  for (i in 2:num_play) {
    for (j in 2:num_play) {    
      
      if (int_m[i,j]=="0"){ #only for the remaining interactions
        if (i < j) {    
          #draw an interaction pair
          r_draw <- runif(1)
          
          #We don't have 1 & 2 yet, because we haven't considered 'avoid' interaction types
          
          #N:N
          if (r_draw < prob_line[3]) {
            #Obtain the interaction types from the probability line, randomize, and assign
            int_name <- unlist(strsplit(strsplit(names(prob_line)[3],"_")[[1]][2],""))
            mij <- sample(int_name,2,replace=FALSE)
            int_m[i,j] <- mij[1]
            int_m[j,i] <- mij[2]
          }
          
          #N:M
          if ((r_draw > prob_line[3]) && (r_draw < prob_line[4])) {
            #Obtain the interaction types from the probability line, randomize, and assign
            int_name <- unlist(strsplit(strsplit(names(prob_line)[4],"_")[[1]][2],""))
            mij <- sample(int_name,2,replace=FALSE)
            int_m[i,j] <- mij[1]
            int_m[j,i] <- mij[2]
          }
          
          #N:I
          if ((r_draw > prob_line[4]) && (r_draw < prob_line[5])) {
            #Obtain the interaction types from the probability line, randomize, and assign
            int_name <- unlist(strsplit(strsplit(names(prob_line)[5],"_")[[1]][2],""))
            mij <- sample(int_name,2,replace=FALSE)
            int_m[i,j] <- mij[1]
            int_m[j,i] <- mij[2]
          }
          
          #I:I
          if ((r_draw > prob_line[5]) && (r_draw < prob_line[6])) {
            #Obtain the interaction types from the probability line, randomize, and assign
            int_name <- unlist(strsplit(strsplit(names(prob_line)[6],"_")[[1]][2],""))
            mij <- sample(int_name,2,replace=FALSE)
            int_m[i,j] <- mij[1]
            int_m[j,i] <- mij[2]
          }
          
        }
      }
    }
  }
  
  #Filling in the diagonal
  #Determine which are active players and which aren't
  for (i in 1:num_play) {
    if (length(which(int_m[i,] == "e")) > 0) {
      int_m[i,i] <- "n"
    } else {
      int_m[i,i] = "i"
    }
  }
  
  #===============================
  #Implement posthoc Rules
  #===============================
  
  #1: Row/Col 1 is the sun
  int_m[1,] <- rep("i",num_play)
  #Things can only eat the sun... they cannot avoid, need, or make it
  #suggestion:
  int_m[which(int_m[,1]!="e"),1]="i"
  
  #2: if player A contains an "m" with player B, player B is "i" with everything
  # except if it is "n" with player A
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
    for (l in 2:(num_play)){
      #To ensure we aren't comparing 3:3, 4:4, etc.
      if (k != l) {
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
  }
  
  #Do we need to run the below code again? I don't think that we do...
  
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
