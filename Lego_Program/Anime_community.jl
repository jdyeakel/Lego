using Distributions
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/anime_sims_func.jl")





#How long does a community take to fill
rep=100;
tend = Array{Int64}(rep);

num_play = 100
init_probs = [
p_n=0.01,
p_a=0.01,
p_m=0.01,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
a_thresh=0.0;
n_thresh=0.0;
for i=1:rep
  int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);
  tout, c_t, c_tp, c_tind = anime_sims_func(int_m,tp_m,tind_m,a_thresh,n_thresh);
  tend[i] = tout;
  print("Rep= ",i)
end
plot(x=tend,Geom.histogram)
