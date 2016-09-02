using Distributions
include("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 30

init_probs = [
p_n=0.02,
p_a=0.2,
p_m=0.1,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m=build_template_degrees(num_play,init_probs)

#Food web network
trophic_m = Array(Int,num_play,num_play)*0;
trophic = find(x->x=='a',int_m);
trophic_m[trophic]=1;
trophic_m
