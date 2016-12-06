using Distributions
# using Gadfly
using RCall
using HDF5
using JLD

@everywhere using Distributions
#@everywhere using Gadfly
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/reppak.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")



S = 400;
reps = 100;
tmax = 1000;

a_thresh = 0;
trophicload = 2;
ppweight = 1/4;
rate_col = 1;

varvec = [0.0001 0.002];
numsp = Array(Array{Int64},length(varvec));
potcol = Array(Array{Int64},length(varvec));
for i=1:length(varvec)
  println("i=",i)
  
  n_thresh = 0.2;

  #Establish community template
  probs = [
  p_n=0.004,
  p_a=0.01,
  p_m= varvec[i], #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]

  #int_m, 
  sprich, 
  rich, 
  conn, 
  #comgen, 
  ext_prim, 
  ext_sec, 
  pot_col,
  num_sp = reppak(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  potcol[i] = pot_col;
  numsp[i] = num_sp;
  
end


# quit()


namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcolS.pdf"
R"""
library(RColorBrewer)
cols = brewer.pal(max(3,length($varvec)),'Set1')
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
pdf($namespace,width=12,height=6);
par(mfrow=c(1,2));
pcol = $(potcol[1]);
plot(seq(1,$tmax)/$(numsp[1][1]),pcol[,1]/$(numsp[1][1]),type='l',xlim=c(0,1),ylim=c(0,0.3),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1])
for (r in 2:$reps) {
  ns = $(numsp[1]);
  lines(seq(1,$tmax)/ns[r],pcol[,r]/ns[r],col=cols[1])
}
legend(x=0.85,y=0.3,legend=$varvec,col=cols,pch=16,title='varvec',cex=1,bty='n');
"""
for i=2:length(varvec)
  R"""
  pcol = $(potcol[i]);
  ns = $(numsp[i]);
  for (r in 1:$reps) {
    lines(seq(1,$tmax)/ns[r],pcol[,r]/ns[r],col=cols[$i])
  }
  """
end 

R"""
nichediffs=c($(diff(potcol[1][2:340,:])),$(diff(potcol[2][2:340,:])));
hist(nichediffs,
breaks=10,xlim=c(0,max(nichediffs)),col=cols[2],xlab='Niche space expansion',main='')
dev.off()
"""




#Without Normalization
R"""
library(RColorBrewer)
cols = brewer.pal(max(3,length($varvec)),'Set1')
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
par(mfrow=c(1,2));
pcol = $(potcol[1]);
plot(seq(1,$tmax),pcol[,1],type='l',xlim=c(0,500),ylim=c(0,100),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1])
for (r in 2:$reps) {
  ns = $(numsp[1]);
  lines(seq(1,$tmax),pcol[,r],col=cols[1])
}
legend(x=0.85,y=0.3,legend=$varvec,col=cols,pch=16,title='varvec',cex=1,bty='n');
"""
for i=2:length(varvec)
  R"""
  pcol = $(potcol[i]);
  ns = $(numsp[i]);
  for (r in 1:$reps) {
    lines(seq(1,$tmax),pcol[,r],col=cols[$i])
  }
  """
end 

R"""
nichediffs=c($(diff(potcol[1][2:340,:])),$(diff(potcol[2][2:340,:])));
hist(nichediffs,
breaks=10,xlim=c(0,max(nichediffs)),col=cols[2],xlab='Niche space expansion',main='')
"""
