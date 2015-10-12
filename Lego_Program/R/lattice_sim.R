lattice_sim <- function(int_m,L){
  source("R/test_compatability.R")
  
  #Lattice size
  size <- (L-1)^2
  
  comm_latt <- list()
  
  #Populate lattice with initial communities
  for (i in 1:size) {
    out <- test_compatability(int_m, N_size=init_size, min_size=10) #Need to add a counter to infer feasibility
    comm_latt[[i]]<-out[[1]]
  }
  
  #Initiate spatial simulations
  
  
  
  
  
  
  
  
  
}