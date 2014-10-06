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
innov.rate <- 5

####################################
#Initial conditions (t = 1)
####################################

#Global properties
niche.space <- numeric(t.term)
sp.richness <- numeric(t.term)
avg.complexity <- numeric(t.term)

#Vector of species
a <- numeric()
a.ab <- numeric()
#List of coproducts
b <- list()
b.ab <- list()
#List of resources
c <- list()
c.ab <- list()

#Define eden species
a[1] <- paste("sp.",rstring(1,str.length),sep="")

#Create vector of initial resource requirements for eden species
c0.size <- round(runif(1,1,10),0)
#Draw random initial resources
c[[1]] <- paste("r.",rstring(c0.size,str.length),sep="")

#Initiate Global Resources + abundances
#The identity of the resource/coproduct/species needs to be paired with numerical abundance
R.id <- c(c[[1]],a[1])
#Draw the abundances of the current resources
R.ab <- sample(seq(1,max.ab),length(R.id)-1)
a1.ab <- 0
R <- data.frame(R.id,c(R.ab,a1.ab),row.names=NULL)
colnames(R) <- c("ID","Abund")


#Establish abundance of Resource use for the eden organism
#i.e. how much of the global resources is the eden organism using?
#This could be a function of an exponential distribution... 
#so that generally resource use is small rather than large
c.ab[[1]] <- numeric(length(c[[1]]))
for (i in 1:length(c[[1]])) {
  res <- which(R$ID == c[[1]][i])
  #Choose whatever is smaller - a random draw or the global amount of each resource
  c.ab[[1]][i] <- min(rexp(1,rate=0.25),R$Abund[res])
}
#Set species abundance
a1.ab <- runif(1)*sum(c.ab[[1]])
a.ab[1] <- a1.ab
#Record species abundance in R
R$Abund[length(c[[1]])+1] <- a1.ab

#Build the coproducts
b[[1]] <- sample(c[[1]],round(runif(1,1,length(c[[1]])-1),0))
#What is the probability that a unique coproduct is formed?
pr.newco <- min(innov.rate * (1/c0.size),1)

#For each coproduct draw a probability of creating a new coproduct
#Create new coproducts in accordance to this probability
draw.newco <- runif(length(b[[1]])) < pr.newco
num.newco <- length(which(draw.newco))
if (num.newco > 0) {
  
  newco.id <- paste("co.",rstring(num.newco,str.length),sep="")
  
  #Update the coproduct list for the eden organism
  b[[1]] <- c(b[[1]],newco.id)

}

#Amount of coproducts PRODUCED must equal the sum(resources) - species abundance
co.prop <- runif(length(b[[1]]))
co.prop <- co.prop/sum(co.prop)
b.ab[[1]] <- co.prop*(sum(c.ab[[1]]) - a1.ab)


#Add 0 amt of new coproducts to the Resource List 
#(just the ID, essentially - ab will be added during next step)
R.id <- c(as.character(R$ID),newco.id)
R.ab <- c(R$Abund,rep(0,length(newco.id)))
R <- data.frame(R.id,R.ab,row.names=NULL)
colnames(R) <- c("ID","Abund")

#Add coproduct abundances to the global resource matrix
#Rebuild Global resource matrix
b.id <- as.numeric(unlist(sapply(b[[1]],function(x){which(x == as.character(R$ID))})))
R$Abund[b.id] <- R$Abund[b.id] + b.ab[[1]]


#Create a Resource-In-Use dataframe
R.inuse <- R
tot.res.use <- numeric(length(R$ID))
for (i in 1:length(R$ID)) {
  res <- as.character(R$ID[i])
  #Find which species is using resource i
  sp.res.use <- which(unlist(lapply(c,function(x){res %in% x})))
  if (length(sp.res.use) > 0) {
    #Obtain the abundances of resource i in use
    sp.res.use.ab <- numeric(length(sp.res.use))
    for (j in 1:length(sp.res.use)) {
      sp.res.use.ab[j] <- c.ab[[sp.res.use[j]]][which(c[[sp.res.use[j]]] == res)]
    }
    tot.res.use[i] <- sum(sp.res.use.ab)
  } else {tot.res.use[i] <- 0}
}
R.inuse$Abund <- tot.res.use

#Note that the coproducts are being PRODUCED. 
#So the number of coproducts in R is the number being produced by the presence of an organism.
#If an organism that produces a coproduct dies - if that coproduct is being used by another organism, 
#it might have to die too



#Calculate Global Properties of the System
niche.space[1] <- length(R$ID)
sp.richness[1] <- length(a)
avg.complexity[1] <- mean(unlist(lapply(c,length)))

trophic.edgelist <- list(t.term)
trophic.edgelist[[1]] <- matrix(0,0,2)
trophic.edgelist.w <- list(t.term)
trophic.edgelist.w[[1]] <- numeric()

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
  new.trophic <- matrix(0,0,2)
  new.trophic.w <- numeric()
  #Initiate mutations ~ new species with small differences in respective b (coproduct) and c (resource) vectors
  if (length(mut.event)>0) {
    
    #The IDs of mut.event determine which species are mutating
    new.sp <- numeric(length(mut.event))
    sp.ab <- numeric(length(mut.event))
    #new.sp.ab <- sample(seq(1,max.ab),length(mut.event))
    mut.res <- list()
    mut.res.ab <- list()
    mut.co <- list()
    mut.co.ab <- list()
    newco.id <- list()
    newco.ab <- list()
    
    for (i in 1:length(mut.event)) {
      #Which species is speciating?
      mutated <- mut.event[i]
      new.sp[i] <- paste("sp.",rstring(i,str.length),sep="")
      #Resources should be similar to resources of mutated species
      mut.res[[i]] <- c[[mutated]]
      mut.res.ab[[i]] <- c.ab[[mutated]]
      
      #Randomly select some deviation ~ we may want to skew these values toward zero (combined -ExpDist, ExpDist)
      mut.res.dev <- runif(length(mut.res[[i]]),-1,1)*max(mut.res.ab[[i]])
      #Apply the deviation
      mut.res.ab[[i]] <- mut.res.ab[[i]] + mut.res.dev
      #Ensure there are no zeros
      elim <- which(mut.res.ab[[i]] <= 0)
      mut.res[[i]] <- mut.res[[i]][-elim]
      mut.res.ab[[i]] <- mut.res.ab[[i]][-elim]
      
      #Add a single new resource? Inventiveness measure
      pr.newres <- min(innov.rate * (1/niche.space[t-1]),1)
      draw.newres <- runif(1) < pr.newco
      
      if (draw.newres) {
        #Ensure that the new resource is not already in it's list
        new.res <- mut.res[[i]][1]
        while (new.res %in% mut.res[[i]] == TRUE) {
          new.res <- as.character(sample(R$ID,1))
        }
        #Define the new resource list for the mutated species
        mut.res[[i]] <- c(mut.res[[i]],new.res)
        
        #Determine the abundance of the new resource
        mut.res.ab[[i]] <- c(mut.res.ab[[i]],min(rexp(1,rate=0.25),
                                                 R$Abund[which(as.character(R$ID) == new.res)] -
                                                   R.inuse$Abund[which(as.character(R$ID) == new.res)]))
      }
      
      #Second, draw a usage between 0 and R-R.use (available resources)
      mut.res.ab[[i]] <- sapply(mut.res.id,function(x){min(rexp(1,rate=0.25),R$Abund[x]-R.inuse$Abund[x])})
      
      #Draw the species abundance
      sp.ab[i] <- runif(1)*sum(mut.res.ab[[i]])
      
      #Record any trophic interactions
      #Identify the prey
      prey <- mut.res[[i]][which(grepl("sp.",mut.res[[i]]))]
      if (length(prey > 0)) {
        trophic <- matrix(0,length(prey),2)
        trophic.w <- numeric()
        for (j in 1:length(prey)) {
          #Record the identities of the trophic interaction
          trophic[j,] <- c(new.sp[i],prey[j])
          #Record the amount consumed
          trophic.w[j] <- mut.res.ab[[i]][which(grepl("sp.",mut.res[[i]]))]
        }
      }
      
      #Determine co-products of new species (a species cannot be a coproduct!)
      co.full <- mut.res[[1]][which(grepl("sp.",mut.res[[1]])==FALSE)]
      mut.co[[i]] <- sample(co.full,round(runif(1,1,length(co.full)-1),0))
      
      #What is the probability that a unique coproduct is formed?
      pr.newco <- min(innov.rate * (1/niche.space[t-1]),1)
      
      #For each coproduct draw a probability of creating a new coproduct
      #Create new coproducts in accordance to this probability
      draw.newco <- runif(length(mut.co[[i]])) < pr.newco
      num.newco <- length(which(draw.newco))
      if (num.newco > 0) {
        
        newco.id[[i]] <- paste("co.",rstring(num.newco,str.length),sep="")
        #Determine the global abundance of the new coproduct
        #newco.ab[[i]] <- sample(seq(1,new.sp.ab[i]/num.newco),length(newco.id[[i]]))
        #newco <- data.frame(newco.id,newco.ab,row.names=NULL)
        
        #Update the coproduct list for the eden organism
        mut.co[[i]] <- c(mut.co[[i]],newco.id[[i]])
        
      }
      
      #Amount of coproducts PRODUCED must equal the sum(resources) - species abundance
      co.prop <- runif(length(mut.co[[i]]))
      co.prop <- co.prop/sum(co.prop)
      mut.co.ab[[i]] <- co.prop*(sum(mut.res.ab[[i]]) - sp.ab[i])
      
      
      #Update trophic interactions IF there is a trophic interaction formed such that length(prey) > 0
      if (length(prey) > 0) {
        if (length(new.trophic[,1]) == 0) {
          new.trophic <- trophic
          new.trophic.w <- trophic.w
        } else {
          new.trophic <- rbind(new.trophic,trophic)
          new.trophic.w <- c(new.trophic.w,trophic.w)
        }
      }
    
    } #End loop over mut.events (i)
    
    #Updating
    
    #Update Species list
    a <- c(a,new.sp)
    a.ab <- c(a.ab,sp.ab)
    c <- c(c,mut.res)
    c.ab <- c(c.ab,mut.res.ab)
    b <- c(b,mut.co)
    b.ab <- c(b.ab,mut.co.ab)
    
    #Add 0 amt of NEW coproducts to the Resource List 
    #(just the ID, essentially - ab will be added during next step)
    R.id <- c(as.character(R$ID),new.sp,unlist(newco.id))
    R.ab <- c(R$Abund,sp.ab,rep(0,length(unlist(newco.id))))
    R <- data.frame(R.id,R.ab,row.names=NULL)
    colnames(R) <- c("ID","Abund")
    
    #Add ALL coproduct abundances to the global resource matrix for each new species i
    #Rebuild Global resource matrix
    for (i in 1:length(mut.event)) {
      b.id <- as.numeric(unlist(sapply(mut.co[[i]],function(x){which(x == as.character(R$ID))})))
      R$Abund[b.id] <- R$Abund[b.id] + mut.co.ab[[i]]
    }    

    #     #Rebuild Global resource matrix IF there are new coproducts AND there are newly evolved species
    #     if ((length(unlist(newco.id)) > 0) && (length(mut.event) > 0)) {
    #       R.id <- c(as.character(R$ID),new.sp,unlist(newco.id))
    #       R.ab <- c(R$Abund,new.sp.ab,unlist(newco.ab))
    #       R <- data.frame(R.id,R.ab,row.names=NULL)
    #       colnames(R) <- c("ID","Abund")
    #     }
    #     #Rebuild Global resource matrix IF there are NO new coproducts AND there are newly evolved species
    #     if ((length(unlist(newco.id)) == 0) && (length(mut.event) > 0)) {
    #       R.id <- c(as.character(R$ID),new.sp)
    #       R.ab <- c(R$Abund,new.sp.ab)
    #       R <- data.frame(R.id,R.ab,row.names=NULL)
    #       colnames(R) <- c("ID","Abund")
    #     }
    
    #Update Resources-In-Use dataframe
    R.inuse.new <- R
    tot.res.use <- numeric(length(R$ID))
    for (i in 1:length(R$ID)) {
      res <- as.character(R$ID[i])
      #Find which species is using resource i
      sp.res.use <- which(unlist(lapply(c,function(x){res %in% x})))
      if (length(sp.res.use) > 0) {
        #Obtain the abundances of resource i in use
        sp.res.use.ab <- numeric(length(sp.res.use))
        for (j in 1:length(sp.res.use)) {
          sp.res.use.ab[j] <- c.ab[[sp.res.use[j]]][which(c[[sp.res.use[j]]] == res)]
        }
        tot.res.use[i] <- sum(sp.res.use.ab)
      } else {tot.res.use[i] <- 0}
    }
    R.inuse.new$Abund <- tot.res.use
    R.inuse <- R.inuse.new
    
  }
  

  #Determine similarities in resource use among species
  #CURRENTLY DOES NOT TAKE ABUNDANCES INTO ACCOUNT
  sim.m <- matrix(0,length(a),length(a))
  for (i in 1:length(a)) {
    for (j in 1:length(a)) {
      sim.m[i,j] <- string.similarity(c[[i]],c.ab[[i]],c[[j]],c.ab[[j]])[2] #[1] = Jaccard; [2] = Cosine Sim Index
    }
  }
  #Eliminate the effect of competition with yourself
  diag(sim.m) <- 0
  comp.pres <- apply(sim.m,2,sum)
  #Need to normalize this measure
  
  #Update system trophic edgelist
  #The consumers are in the left column and the prey are in the right column
  #If there has NOT been an additional trophic interaction formed
  if (length(new.trophic) == 0) {
    trophic.edgelist[[t]] <- trophic.edgelist[[t-1]]
    trophic.edgelist.w[[t]] <- trophic.edgelist.w[[t-1]]
  } else {
    #If this is the first trophic interaction, AND there has been an additional trophic interaction formed
    if ((length(trophic.edgelist[[t-1]][,1]) == 0) && (length(new.trophic) > 0)) {
      trophic.edgelist[[t]] <- new.trophic
      trophic.edgelist.w[[t]] <- new.trophic.w
    } 
    #If this is not the first trophic interaction, and there has been an additional trophic interaction formed
    if ((length(trophic.edgelist[[t-1]][,1]) > 0) && (length(new.trophic) > 0)) {
      trophic.edgelist[[t]] <- rbind(trophic.edgelist[[t-1]],new.trophic)
      trophic.edgelist.w[[t]] <- c(trophic.edgelist[[t-1]],new.trophic.w)
    }
  }
  
  #What are the predation pressures for each species?
  pred.pres <- numeric(length(a))
  for (j in 1:length(a)) {
    pred.pres[j] <- max(0,trophic.edgelist.w[[t]][which(trophic.edgelist[[t]][,2] == a[j])]) / R$Abund[which(as.character(R$ID) == a[j])]
  }
  
  #Set dynamic rates
  #Extinction rate :: should increase with competition load and predation load
  ext.background <- 0.01
  ext.rate <- ext.background + 0.5*comp.pres + 0.5*pred.pres
  draw.ext <- runif(length(a),0,1) < ext.rate
  extinct <- which(draw.ext)
  
  #Induce extinctions
  #Modify R.inuse matrix
  #Modify trophic interactions + trophic interaction weights
  #Determine whether coproducts of extinct species are unique
  #If coproducts are unique, eliminate them from the Global Resources
  
  
  
  #Delete trophic interactions involving the extinct species
  pred.extinct <- which(trophic.edgelist[[t]][,1] == a[extinct])
  prey.extinct <- which(trophic.edgelist[[t]][,2] == a[extinct]) #this also records the ID of the predators that are consuming the extinct prey
  trophic.edgelist[[t]] <- trophic.edgelist[[t]][-c(pred.extinct,prey.extinct),]
  trophic.edgelist.w[[t]] <- trophic.edgelist.w[[t]][-c(pred.extinct,prey.extinct)]
  
  #Modify Resources-in-use by extinct species
  modify.res.id <- unlist(c[extinct])
  modify.res.ab <- unlist(c.ab[extinct])
  
  for (j in 1:length(prey.extinct)) {
    pred.id <- prey.extinct[j]
    prey.id <- 
    c[[pred.id]] <- c[[pred.id]][-which(c[[pred.id]] == )]
    c.ab[[j]] 
  }
  
  
  
  #Delete species from the a vector
  a <- a[-extinct]
  
  
  
}




























