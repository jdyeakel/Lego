#########################
# Niche Evolution Model #
#########################

source("rstring.r")
str_length <- 10

#Initial conditions

#Vector of species
a <- numeric()
#List of coproducts
b <- list()
#List of resources
c <- list()

#Define eden species
a[1] <- paste("sp_",rstring(1,str_length),sep="")

#Create vector of initial resource requirements for eden species
c0_size <- round(runif(1,0,10),0)
#Draw random initial resources
c[[1]] <- rstring(c0_size,str_length)

#Build the coproducts
#What is the probability that a unique coproduct is formed?
pr_newco <- 1/length(c0)
b[[1]] <- sample(c[[1]],round(runif(1,0,length(c[[1]])),0))
#For each coproduct draw a probability of creating a new coproduct
#Create new coproducts in accordance to this probability
draw_newco <- runif(length(b[[1]])) < pr_newco
newco <- rstring(length(which(draw_newco)),str_length)
b[[1]] <- c(b[[1]],newco)