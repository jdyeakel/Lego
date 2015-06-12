rm(list=c(ls()))


num_play <- 50

#Interaction probabilities
#Must add to 1
pr_a <- 1/5
pr_n <- 1/5
pr_i <- 1/5
pr_m <- 1/5
pr_e <- 1/5

#Pairwise interaction probabilities

# pr_aa <- pr_a*(pr_a/(pr_a + pr_i))
# pr_ai <- pr_a*(pr_i/(pr_a + pr_i))
# 
# pr_nn <- pr_n*(pr_n/(pr_i + pr_n + pr_e + pr_m))
# pr_ni <- pr_n*(pr_i/(pr_i + pr_n + pr_e + pr_m))
# pr_nm <- pr_n*(pr_m/(pr_i + pr_n + pr_e + pr_m))
# pr_ne <- pr_n*(pr_e/(pr_i + pr_n + pr_e + pr_m))
# 
# pr_ia <- pr_i*(pr_a/(pr_a + pr_n + pr_i + pr_e))
# pr_in <- pr_i*(pr_n/(pr_a + pr_n + pr_i + pr_e))
# pr_ii <- pr_i*(pr_i/(pr_a + pr_n + pr_i + pr_e))
# pr_ie <- pr_i*(pr_e/(pr_a + pr_n + pr_i + pr_e))
# 
# pr_mn <- pr_m*pr_n
# 
# pr_ei <- pr_e*(pr_i/(pr_i + pr_n + pr_e))
# pr_en <- pr_e*(pr_n/(pr_i + pr_n + pr_e))
# pr_ee <- pr_e*(pr_e/(pr_i + pr_n + pr_e))

pw_prob <- c(
  pr_ne = 0.025,
  pr_nn = 0.025,
  pr_ni = 0.05,
  pr_nm = 0.05,
  pr_ia = 0.05,
  pr_ie = 0.2,
  pr_ii = 0.5,
  pr_aa = 0.05,
  pr_ee = 0.05
)
#make sure this vector sums to one
pw_prob <- pw_prob / sum(pw_prob)

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


#for 'need' coexistence: if vector of need for spA != player vector, spA is locally extinct


