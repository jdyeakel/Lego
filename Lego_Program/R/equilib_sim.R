equilib_sim <- function(int_m,init_size,reps){
  source("R/test_compatability.R")
  
  comm_comp <- list()
  comm_dim <- numeric(reps)
  #Community trophic degree distribution
  comm_dd <- list()
  
  for (r in 1:reps) {
    #Build sub-communities from master template
    out<-test_compatability(int_m, N_size=init_size, min_size=10) #Need to add a counter to infer feasibility
    sub_m<-out[[1]]
    #feasibility=out[[2]] #iterations to find a feasible subcommunity
    comm_comp[[r]] <- rownames(sub_m)
    comm_dim[r] <- dim(sub_m)[1]
    
    #trophic degree distributions
    trophic_dd <- numeric(dim(sub_m)[1])
    for (i in 2:dim(sub_m)[1]) {
      trophic_dd[i] <- length(which(sub_m[,i] == "e"))
    }
    comm_dd[[r]] <- trophic_dd
  }
  
  
  
  #Compare the similarity of assembled communities
  #One way is to just compare the presence/absence of each row/column
  comm_sim <- matrix(0,reps,reps)
  for (i in 1:reps) {
    for (j in 1:reps) {
      if (i <= j) {
        sp_i <- comm_comp[[i]]
        sp_j <- comm_comp[[j]]
        
        comm_sim[i,j] <- length(intersect(sp_i,sp_j))/length(union(sp_i,sp_j))
        
        
      }
    }
    
  }
  
  Rout <- list()
  Rout[[1]] <- comm_sim
  Rout[[2]] <- comm_dd
  return(Rout)
  
}