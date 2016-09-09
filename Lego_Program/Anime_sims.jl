using Distributions
using Gadfly
include("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 1000

init_probs = [
p_n=0.05,
p_a=0.02,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m, t_m = build_template_degrees(num_play,init_probs);
writedlm("$(homedir())/Dropbox/PostDoc/2014_Lego/data_template/fw.csv",t_m);
