using Distributions
using Gadfly
include("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 10000

init_probs = [
p_n=0.01,
p_a=0.01,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m, sp_m, t_m, tp_m = build_template_degrees(num_play,init_probs);
writedlm("$(homedir())/Dropbox/PostDoc/2014_Lego/data_template/fw.csv",t_m);


#Order of operations
#1) Establish thresholds
#2) Choose root species (must be primary producer)
#3) Choose next species...
  #3a) determine and test 'a' and 'n' thresholds
  #3b) assess similarity to determine exclusion
  #3c) Pass/Fail
#4) Cumulatively add species over time, where inclusion is determined by 1) thresholds, 2) exclusion


#NOTE: NEED a matrix of JUST SPECIES interactions

#primary species
prim_prod = find(x->x=='a',int_m[:,1]);
draw = rand(prim_prod);
seed = int_m[draw,:];
seedrev = int_m[:,draw];
seedmake = find(x->x=='m',seed);
init_play = vcat(seed,int_m[seedmake,:])

tmax = 100;

for t=1:tmax

  #Current species list
  if t==1
    slist = seed;
    srevlist = seedrev;
  end

  #draw a random species
  #continue to draw until you find a trophically-able species
  #JUST GO THROUGH SPECIES!
  contrun=true;
  samplespp = sample(S,collect(1:S),replace=false);
  while contrun == true
    #Draw a new species from the pool
    #List of things that they make
    #Do they consume species or things they make that are present?
