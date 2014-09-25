#########################
# Niche Evolution Model #
#########################

source("rstring.r")
source("insertRow.r")

str_length <- 10
max_ab <- 1000

#Initial conditions

#Vector of species
a <- numeric()
#List of coproducts
b <- list()
#List of resources
c <- list()
c_ab <- list()

#Define eden species
a[1] <- paste("sp_",rstring(1,str_length),sep="")

#Create vector of initial resource requirements for eden species
c0_size <- round(runif(1,0,10),0)
#Draw random initial resources
c[[1]] <- rstring(c0_size,str_length)

#Initiate Global Resources + abundances
#The identity of the resource/coproduct/species needs to be paired with numerical abundance
Rid <- c(c[[1]],a[1])
Rab <- sample(seq(1,max_ab),length(Rid))
R <- data.frame(Rid,Rab,row.names=NULL)
colnames(R) <- c("ID","Abund")

#Build the coproducts
#What is the probability that a unique coproduct is formed?
pr_newco <- 1/length(c0)
b[[1]] <- sample(c[[1]],round(runif(1,0,length(c[[1]])),0))
#For each coproduct draw a probability of creating a new coproduct
#Create new coproducts in accordance to this probability
draw_newco <- runif(length(b[[1]])) < pr_newco
num_newco <- length(which(draw_newco))
if (num_newco > 0) {
  
  newco_id <- rstring(num_newco,str_length)
  newco_ab <- sample(seq(1,max_ab),length(newco_id))
  #newco <- data.frame(newco_id,newco_ab,row.names=NULL)
  
  #Update the coproduct list for the eden organism
  b[[1]] <- c(b[[1]],newco_id)
  
  #Rebuild Global resource matrix
  Rid <- c(as.character(R$ID),newco_id)
  Rab <- c(R$Abund,newco_ab)
  R <- data.frame(Rid,Rab,row.names=NULL)
  colnames(R) <- c("ID","Abund")

}

#Establish Resource use for the eden organism
c_ab[[1]] <- 






