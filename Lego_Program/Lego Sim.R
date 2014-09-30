#########################
# Niche Evolution Model #
#########################

source("rstring.r")
source("insertRow.r")

str.length <- 10
max.ab <- 1000
t.term <- 100

#Initial conditions (t = 1)

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
pr.newco <- 1/length(c0)

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


#Calculate Global Properties of the System
niche.space[1] <- length(R$ID)
sp.richness[1] <- length(a)
avg.complexity[1] <- mean(unlist(lapply(c,length)))














