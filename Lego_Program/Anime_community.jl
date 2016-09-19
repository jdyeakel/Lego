using Distributions
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/anime_sims_func.jl")

num_play = 200

init_probs = [
p_n=0.01,
p_a=0.01,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]


#How long does a community take to fill
rep=100;
tend = Array{Int64}(rep);
a_thresh=0.0;
n_thresh=0.2;
for i=1:rep
  int_m, tout, c_t = anime_sims_func(num_play,init_probs,a_thresh,n_thresh);
  tend[i] = tout;
  # print("i= ",i)
end
plot(x=tend,Geom.histogram)
