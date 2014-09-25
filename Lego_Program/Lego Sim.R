#########################
# Niche Evolution Model #
#########################

source("rstring.r")
str_length <- 10

#Initial conditions

#Create vector of initial resource requirements for eden species
c0_size <- round(runif(1,0,10),0)
#Draw random initial resources
c0 <- rstring(c0_size,str_length)
#What is the probability that a unique coproduct is formed?
pr_newco <- 1/length(c0)
