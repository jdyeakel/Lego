@everywhere using Distributions
#@everywhere using Gadfly
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/reppak.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")



num_play = 500;
reps = 20;
tmax = 500;

a_thresh = 0;
trophicload = 2;
ppweight = 1/4;
rate_col = 1;

varvec = [0.1 0.2 0.3];
numsp = zeros(varvec);
potcol = Array(Array{Int64},length(varvec));
for i=1:length(varvec)
  println("i=",i)
  
  n_thresh = varvec[i];

  #Establish community template
  probs = [
  p_n=0.004,
  p_a=0.01,
  p_m= 0.002, #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]

  #int_m, 
  sprich, 
  rich, 
  conn, 
  comgen, 
  ext_prim, 
  ext_sec, 
  pot_col = reppak(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  numsp[i] = length(find(x->x=='n',diag(int_m)));
  potcol[i] = pot_col;
  
end


namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol.pdf"
R"""
library(RColorBrewer)
cols = brewer.pal(max(3,length($varvec)),'Set1')
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
pdf($namespace,width=8,height=6);
pcol = $(potcol[1]);
plot(seq(1,$tmax)/$(numsp[1]),pcol[,1]/$(numsp[1]),type='l',xlim=c(0,1),ylim=c(0,0.3),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1])
for (r in 2:$reps) {
  lines(seq(1,$tmax)/$(numsp[1]),pcol[,r]/$(numsp[1]),col=cols[1])
}
legend(x=0.9,y=0.3,legend=$varvec,col=cols,pch=16,title='varvec',cex=1,bty='n');
"""
for i=2:length(varvec)
  R"""
  pcol = $(potcol[i]);
  for (r in 1:$reps) {
    lines(seq(1,$tmax)/$(numsp[i]),pcol[,r]/$(numsp[i]),col=cols[$i])
  }
  """
end 
R"dev.off()"
