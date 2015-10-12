
rm(list=c(ls()))
#setwd("/Users/piresmm/git/Lego/Lego_Program/")
library(igraph)
library(plotrix)
library(RColorBrewer)

source("R/build_template_degrees.R")
source("R/test_compatability.R")

#sequence <- seq(10,2000,100)
#prop_active <- numeric(length(sequence))


tic <- 0

i<-50 #template size
#for (i in sequence) {

tic <- tic + 1

num_play <- i

#Defining probabilities of each type
#Basic probs could also be based on num_play. e.g., We should expect p.n*num_play n's per column/row

p.n=0.02 
p.e=0.05
p.m=0.1
p.a=0
p.i= 1 - (sum(p.n,p.e,p.m,p.a))#Ignore with 1 - pr(sum(other))

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
#===============================

#Visualize matrix:
plot_matrix<-function(int_m, num_play){
xx=matrix(as.numeric(as.factor(int_m)),c(num_play,num_play))
par(mar=c(1,1,1,4))
int_types=levels(as.factor(int_m))
color2D.matplot(xx,extremes=c(1:length(int_types)), border="white", axes=FALSE, xlab="", ylab="",main="")
legend(x=num_play+1,y=num_play,legend=int_types,pch=22,pt.bg=c(1:5),xpd=TRUE, bty="n")
}
plot_matrix(int_m, num_play)

#Extract food web
pal <- brewer.pal(3,"Set1")
fw <- (int_m == "e")*1
fw <- fw[-which(player_id==0),-which(player_id==0)]
fw_g <- graph.adjacency(fw)
basal_pos <- 1
trophic <- sapply(1:vcount(fw_g),function(x){shortest.paths(fw_g,v = basal_pos, to = x)})
trophic[which(trophic==Inf)] <- 0
coords <- cbind(runif(vcount(fw_g)),trophic); coords[basal_pos,] <- c(0.5,trophic[basal_pos])
par(mar=c(1,1,1,1))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.0,
     main=ecount(fw_g)/num_play^2,vertex.label=NA,
     vertex.color=pal[2])

###############################

#Compatibility of subsampled community

#for 'need' coexistence: if vector of need for spA != player vector, spA is locally extinct

source("R/test_compatability.R")

N_size=30 #sample richness
min_size=10 #minimum number of elements

#Lattice edge size
#L determines the size for the internal lattice.
#The bordering rows and columns
L <- 5 
size <- (L+2)^2
#Nearest neighbor function for Torus
source("R/ipbc.R")

#Populate the metacommunity with local interaction webs
meta_com <- as.list(rep(0,size))
for (i in 1:size) {
  local <- test_compatability(int_m, N_size, min_size)
  meta_com[[i]] <- local
}

#Interactions...
#===============================
#Topology - Template vs subsample
#==================================

#Generating template
num_play=200
int_m <- build_template(num_play,pw_prob, 0.8)
plot_matrix(int_m, num_play)

#Generating subsample
out<-test_compatability(int_m, N_size=100, min_size=30) #Need to add a counter to infer feasibility
sub_m<-out[[1]]
feasibility=out[[2]] #iterations to find a feasible subcommunity
plot_matrix(sub_m, dim(sub_m)[1])


#.........................
#proportion of active players
#.........................

a_players=sum(diag(int_m)=="n")
in_players=sum(diag(int_m)=="i")

prop.active_temp<-a_players/num_play

a_players_sub=sum(diag(sub_m)=="n")
in_players_sub=sum(diag(sub_m)=="i")

prop.active_sub<-a_players_sub/nrow(sub_m)

#................................
#proportion of interaction types
#...............................
prop.interactios_temp<-c(E=sum(int_m=="e"),
                    N=sum(int_m=="n")-a_players,
                    I=sum(int_m=="i")-in_players,
                    M=sum(int_m=="m"),
                    A=sum(int_m=="a")) /(num_play*(num_play-1))
                  
prop.interactios_sub<-c(E=sum(sub_m=="e"),
                        N=sum(sub_m=="n")-a_players_sub,
                        I=sum(sub_m=="i")-in_players_sub,
                        M=sum(sub_m=="m"),
                        A=sum(sub_m=="a")) /(num_play*(num_play-1))

#...................
#Degree distribution
#Separate in and out degree?
#...................
#template
a_play<-which(diag(int_m)=="n") #active players
species_m<-int_m[a_play,a_play] #matrix with species only
meank<-mean(apply(species_m=="e",1,sum))

#subset
a_play<-which(diag(sub_m)=="n")#active players
species_m<-sub_m[a_play,a_play]#matrix with species only
meank<-mean(apply(species_m=="e",1,sum))

#plotting
par(mfrow=c(2,1),mar=c(3,3,3,4)) 
hist(apply(species_m=="e",1,sum), main=paste("k = ",round(meank,1)), xlab="Degree", col="darkgrey")
hist(apply(species_m=="e",1,sum), main=paste("k = ",round(meank,1)), xlab="Degree", col="darkgrey")















