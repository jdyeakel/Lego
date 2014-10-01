#########################
# Niche Evolution Model #
#########################

rm(list=c(ls()))


source("rstring.r")
source("insertRow.r")
source("string.similarity.r")

str.length <- 10
max.ab <- 1000
t.term <- 100

####################################
#Initial conditions (t = 1)
####################################

#Global properties
niche.space <- numeric(t.term)
sp.richness <- numeric(t.term)
avg.complexity <- numeric(t.term)

#Vector of species
a <- numeric()
#List of coproducts
b <- list()
#List of resources
c <- list()
c.ab <- list()

#Define eden species
a[1] <- paste("sp.",rstring(1,str.length),sep="")

#Create vector of initial resource requirements for eden species
c0.size <- round(runif(1,0,10),0)
#Draw random initial resources
c[[1]] <- rstring(c0.size,str.length)

#Initiate Global Resources + abundances
#The identity of the resource/coproduct/species needs to be paired with numerical abundance
R.id <- c(c[[1]],a[1])
#Draw the abundances of the current resources
R.ab <- sample(seq(1,max.ab),length(R.id))
R <- data.frame(R.id,R.ab,row.names=NULL)
colnames(R) <- c("ID","Abund")

#Build the coproducts
b[[1]] <- sample(c[[1]],round(runif(1,1,length(c[[1]])),0))
#What is the probability that a unique coproduct is formed?
pr.newco <- 1/c0.size

#For each coproduct draw a probability of creating a new coproduct
#Create new coproducts in accordance to this probability
draw.newco <- runif(length(b[[1]])) < pr.newco
num.newco <- length(which(draw.newco))
if (num.newco > 0) {
  
  newco.id <- rstring(num.newco,str.length)
  #Determine the global abundance of the new coproduct
  newco.ab <- sample(seq(1,max.ab),length(newco.id))
  #newco <- data.frame(newco.id,newco.ab,row.names=NULL)
  
  #Update the coproduct list for the eden organism
  b[[1]] <- c(b[[1]],newco.id)
  
  #Rebuild Global resource matrix
  R.id <- c(as.character(R$ID),newco.id)
  R.ab <- c(R$Abund,newco.ab)
  R <- data.frame(R.id,R.ab,row.names=NULL)
  colnames(R) <- c("ID","Abund")

}

#Establish Resource use for the eden organism
#i.e. how much of the global resources is the eden organism using?
#This could be a function of an exponential distribution... so that generally resource use is small rather than large
c.ab[[1]] <- numeric(length(c[[1]]))
for (i in 1:length(c[[1]])) {
  res <- which(R$ID == c[[1]][i])
  #Choose whatever is smaller - a random draw or the global amount of each resource
  c.ab[[1]][i] <- min(rexp(1,rate=0.25),R$Abund[res])
}

#Create a Resource-In-Use dataframe
R.inuse <- R
tot.res.use <- numeric(length(R$ID))
for (i in 1:length(R$ID)) {
  res <- as.character(R$ID[i])
  #Find which species is using resource i
  sp.res.use <- unlist(lapply(c,function(x){which(x==res)}))
  if (length(sp.res.use) > 0) {
    #Obtain the abundances of resource i in use
    sp.res.use.ab <- numeric(length(a))
    for (j in 1:length(a)) {
      sp.res.use.ab[j] <- c.ab[[j]][sp.res.use[j]]
    }
    tot.res.use[i] <- sum(sp.res.use.ab)
  } else {tot.res.use[i] <- 0}
}
R.inuse$Abund <- tot.res.use


#Calculate Global Properties of the System
niche.space[1] <- length(R$ID)
sp.richness[1] <- length(a)
avg.complexity[1] <- mean(unlist(lapply(c,length)))





####################################
#Forward Iterations
####################################


#Dynamic properties of the system




#Begin time iterations
for (t in 2:t.term) {
  
  #First, allow mutations for each species
  #Mutation rate
  m.rate <- 0.01
  m.draw <- runif(length(a),0,1) < m.rate
  mut.event <- which(m.draw)
  
  #Initiate mutations ~ new species with small differences in respective b (coproduct) and c (resource) vectors
  if (length(mut.event)>0) {
    
    #The IDs of mut.event determine which species are mutating
    new.sp <- numeric(length(mut.event))
    mut.res <- list(length(mut.event))
    mut.res.ab <- list(length(mut.event))
    mut.co <- list(length(mut.event))
    newco.id <- list(length(mut.event))
    newco.ab <- list(length(mut.event))
    
    for (i in 1:length(mut.event)) {
      #Which species is speciating?
      mutated <- mut.event[i]
      new.sp[i] <- paste("sp.",rstring(i,str.length),sep="")
      #Resources should be similar to resources of mutated species
      old.res <- c[[mutated]]
      
      #Build new resource list. Should be similar, but not the same.
      #For now, we can just delete a random resource and replace it with a random draw from the R list
      mut.res[[i]] <- old.res[-sample(seq(1,length(old.res)),1)]
      
      #Ensure that the new resource is not already in it's list
      new.res <- mut.res[[i]][1]
      while (new.res %in% mut.res[[i]] == TRUE) {
        new.res <- as.character(sample(R$ID,1))
      }
      #Define the new resource list for the mutated species
      mut.res[[i]] <- c(mut.res[[i]],new.res)
      
      #Draw resource use for mutated species... must be between [0 : R-R.inuse]
      #First get the location of the resource in R
      mut.res.id <- as.numeric(sapply(mut.res[[i]],function(x){which(x==as.character(R$ID))}))
      #Second, draw a usage between 0 and R-R.use (available resources)
      mut.res.ab[[i]] <- sapply(mut.res.id,function(x){min(rexp(1,rate=0.25),R$Abund[x]-R.inuse$Abund[x])})
    
      #Determine co-products of new species
      #Build the coproducts
      mut.co[[i]] <- sample(as.character(R$ID),round(runif(1,1,length(as.character(R$ID))),0))
      while(any(grepl("sp.",mut.co[[i]]))) {
        mut.co[[i]] <- sample(as.character(R$ID),round(runif(1,1,length(as.character(R$ID))),0))
      }
      #What is the probability that a unique coproduct is formed?
      pr.newco <- 1/niche.space[t-1]
      
      #For each coproduct draw a probability of creating a new coproduct
      #Create new coproducts in accordance to this probability
      draw.newco <- runif(length(mut.co[[i]])) < pr.newco
      num.newco <- length(which(draw.newco))
      if (num.newco > 0) {
        
        newco.id[[i]] <- rstring(num.newco,str.length)
        #Determine the global abundance of the new coproduct
        newco.ab[[i]] <- sample(seq(1,max.ab),length(newco.id))
        #newco <- data.frame(newco.id,newco.ab,row.names=NULL)
        
        #Update the coproduct list for the eden organism
        mut.co[[i]] <- c(mut.co[[i]],newco.id[[i]])
        
      }
    
    } #End loop over mut.events
    
    #Update Species list
    a <- c(a,new.sp)
    c <- c(c,mut.res)
    c.ab <- c(c.ab,mut.res.ab)
    b <- c(b,mut.co)
    
    #Rebuild Global resource matrix IF there are new coproducts AND there are newly evolved species
    if ((length(unlist(newco.id)) > 0) && (length(mut.event) > 0)) {
      R.id <- c(as.character(R$ID),new.sp,unlist(newco.id))
      R.ab <- c(R$Abund,sample(seq(1,max.ab),length(mut.event)),unlist(newco.ab))
      R <- data.frame(R.id,R.ab,row.names=NULL)
      colnames(R) <- c("ID","Abund")
    }
    #Rebuild Global resource matrix IF there are NO new coproducts AND there are newly evolved species
    if ((length(unlist(newco.id)) == 0) && (length(mut.event) > 0)) {
      R.id <- c(as.character(R$ID),new.sp)
      R.ab <- c(R$Abund,sample(seq(1,max.ab),length(mut.event)))
      R <- data.frame(R.id,R.ab,row.names=NULL)
      colnames(R) <- c("ID","Abund")
    }
    
    #Update Resources-In-Use dataframe
    
    
  }
  
  
  
  
  #Determine similarities in resource use among species
  #CURRENTLY DOES NOT TAKE ABUNDANCES INTO ACCOUNT
  sim.m <- matrix(0,length(a),length(a))
  for (i in 1:length(a)) {
    for (j in 1:length(a)) {
      sim.m[i,j] <- string.similarity(c[[i]],c[[j]])
    }
  }
  
  #Set dynamic rates
  
  
  #Extinction rate
  #For each species, what is the degree of competition?
  
  #For each species, what is the degree of predation?
  
  ext.rate <- 
  
  
  
  
  
  #For each species, determine whether a mutation occurs
  
  
  
}




























