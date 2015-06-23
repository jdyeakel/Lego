
rm(list=c(ls()))
#setwd("/Users/piresmm/git/Lego/Lego_Program/")
library(igraph)
library(plotrix)
library(RColorBrewer)

source("R/build_template.R")
source("R/test_compatability.R")

#sequence <- seq(10,2000,100)
#prop_active <- numeric(length(sequence))


tic <- 0

i<-50 #template size
#for (i in sequence) {

tic <- tic + 1

num_play <- i

#   pw_prob <- c(
#     pr_ne = 0.025,
#     pr_nn = 0.025,
#     pr_ni = 0.05,
#     pr_nm = 0.005,
#     pr_ia = 0.05,
#     pr_ie = 0.2,
#     pr_ii = 0.5,
#     pr_aa = 0.05,
#     pr_ee = 0.05
#   )

#Defining probabilities of each type
#Basic probs could also be based on num_play. e.g., We should expect p.n*num_play n's per column/row

p.n=0.02 
p.e=0.1
p.i=0.72
p.m=0.1
p.a=0

#Normalization [0,1]
S_prob=sum(c(p.n,p.e,p.i,p.m,p.a))
p.n=p.n/S_prob
p.e=p.e/S_prob
p.i=p.i/S_prob
p.m=p.m/S_prob
p.a=p.a/S_prob

#Defining paiwise probabilities 
pw_prob <- c(
  pr_ne = p.n*(p.e/(p.e+p.n+p.i+p.m)),
  pr_nn = p.n*(p.n/(p.e+p.n+p.i+p.m)),
  pr_ni = p.n*(p.i/(p.e+p.n+p.i+p.m)),
  pr_nm = p.n*(p.m/(p.e+p.n+p.i+p.m)),
  pr_ia = p.i*(p.a/(p.e+p.a+p.n+p.i)),
  pr_ie = p.i*(p.e/(p.e+p.a+p.n+p.i)),
  pr_ii = p.i*(p.i/(p.e+p.a+p.n+p.i)),
  pr_aa = p.a*(p.a/(p.a+p.i)),
  pr_ee = p.e*(p.e/(p.i+p.n+p.e))
)


#make sure this vector sums to one
pw_prob <- pw_prob / sum(pw_prob)

#Build the interaction template
int_m <- build_template(num_play,pw_prob, 0.8)

#Some funky things
# Passive players can 'make' things... this might be okay (chemical rxns)

#Proportion of active players
player_id <- apply(int_m,1,function(x){length(which(x == "e")) > 0})*1
#prop_active[tic] <- sum(player_id)/num_play

#fw_connectance[tic] <- 


#}

#plot(sequence,prop_active,xlab="size",ylab="proportion active players",ylim=c(0,1))


#VISUALIZATION OF THE TEMPLATE
###############################

#Visualize matrix:
xx=matrix(as.numeric(as.factor(int_m)),c(num_play,num_play))
par(mar=c(1,1,1,4))
color2D.matplot(xx,extremes=c(1:5), border="white", axes=FALSE, xlab="", ylab="",main="")
legend(i,i,legend=levels(as.factor(int_m)),pch=22,pt.bg=c(1:5),xpd=TRUE, bty="n")

#Extract food web
pal <- brewer.pal(3,"Set1")
fw <- (int_m == "e")*1
fw <- fw[-which(player_id==0),-which(player_id==0)]
fw_g <- graph.adjacency(fw)
basal_pos <- 1
trophic <- sapply(1:vcount(fw_g),function(x){shortest.paths(fw_g,v = basal_pos, to = x)})
trophic[which(trophic==Inf)] <- 0
coords <- cbind(runif(num_play),trophic); coords[basal_pos,] <- c(0.5,trophic[basal_pos])
par(mar=c(1,1,1,1))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.4,
     main=ecount(fw_g)/num_play^2,vertex.label=NA,
     vertex.color=pal[player_id+1])

###############################

#Compatibility of subsampled community

#for 'need' coexistence: if vector of need for spA != player vector, spA is locally extinct

source("R/test_compatability.R")

R<-1 #counter
min_size=10 #minimum number of elements

#iterate until finding a viable community with at least min_size elements
while (R<min_size){ 
  
  num_play<-dim(int_m)[1] #template size
  N_size=30 #sample richness
  r_sample <- c(1,sort(sample(2:num_play,N_size))) #sun + sampled elements 
  
  #Testing
  local=test_compatability (int_m,r_sample)
  R=dim(local)[1] #community size
}















