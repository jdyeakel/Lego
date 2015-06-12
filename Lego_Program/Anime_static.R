rm(list=c(ls()))

library(igraph)
library(plotrix)
library(RColorBrewer)

source("R/build_template.R")

num_play <- 100

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

#Build the interaction template
int_m <- build_template(num_play,pw_prob)

#Some funky things
# Passive players can 'make' things... this might be okay (chemical rxns)

#Proportion of active players
sum(apply(int_m,1,function(x){length(which(x == "e")) > 0})*1)/num_play


#Visualize matrix:
xx=matrix(as.numeric(as.factor(int_m)),c(num_play,num_play))
color2D.matplot(xx,extremes=c(1:5), border="white", axes=FALSE, xlab="", ylab="",main="")

#Extract food web
fw <- (int_m == "e")*1
fw_g <- graph.adjacency(fw)
basal_pos <- 1
trophic <- sapply(1:num_play,function(x){shortest.paths(fw_g,v = basal_pos, to = x)})
trophic[which(trophic==Inf)] <- 0
coords <- cbind(runif(num_play),trophic); coords[basal_pos,] <- c(0.5,trophic[basal_pos])
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.04,main=ecount(fw_g)/num_play^2,vertex.label=NA)

#for 'need' coexistence: if vector of need for spA != player vector, spA is locally extinct


