using Distributions
#using Gadfly
using RCall
using HDF5
using JLD

include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")


S = 400;
ppweight=1/4;

probs = [
p_n=0.004,
p_a=0.01,
p_m=0.003,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_species(S,probs,ppweight);
length(find(x->x=='n',diag(int_m)))

num_play = 500;
int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_degrees(num_play,probs,ppweight);
