using Distributions
include("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 10

init_probs = [
p_n=0.02,
p_a=0.05,
p_m=0.1,
p_e=0.0,
p_i= 1 - sum([p_n,p_e,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

build_template_degrees(num_play,probs)
