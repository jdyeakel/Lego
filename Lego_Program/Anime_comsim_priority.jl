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
reps = 5;
tmax = 500;

a_thresh = 0;
n_thresh = 0.2;
trophicload = 2;
ppweight = 1/4;
rate_col = 1;

#Establish community template
probs = [
p_n=0.004,
p_a=0.01,
p_m= 0.00002, #needvec[i]/num_play, #1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m, sprich, rich, conn, comgen, ext_prim, ext_sec, pot_col = reppak(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

numsp = length(find(x->x=='n',diag(int_m)));
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol.pdf"
R"""
pdf($namespace,width=8,height=6);
pcol = $pot_col;
plot(seq(1,$tmax)/$numsp,pcol[,1]/$numsp,type='l',xlim=c(0,1),ylim=c(0,0.2),xlab='Proportion filled',ylab='Proportion potential colonizers')
for (r in 2:$reps) {
  lines(seq(1,$tmax)/$numsp,pcol[,r]/$numsp)
}
dev.off()
"""
