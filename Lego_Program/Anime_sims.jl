using Distributions
using Gadfly
include("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 10

init_probs = [
p_n=0.01,
p_a=0.1,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);
#writedlm("$(homedir())/Dropbox/PostDoc/2014_Lego/data_template/fw.csv",t_m);


#Order of operations
#1) Establish thresholds
#2) Choose root species (must be primary producer)
#3) Choose next species...
  #3a) determine and test 'a' and 'n' thresholds
  #3b) assess similarity to determine exclusion
  #3c) Pass/Fail
#4) Cumulatively add species over time, where inclusion is determined by 1) thresholds, 2) exclusion


#Setting thresholds
n_thresh = 0.0;
a_thresh = 0.0;


#NOTE: NEED a matrix of JUST SPECIES interactions
Slist = find(x->x=='n',diag(int_m));

#primary species
prim_prod = find(x->x=='a',int_m[:,1]);
id = rand(prim_prod);
#interactions of the seed species
seed = int_m[id,:];
#reverse interactions of the seed species
seedrev = int_m[:,id];
#What things does the seed species make?
idm = find(x->x=='m',seed);
#What things does the seed species need?
idn = find(x->x=='n',seed);

#Interactions of idseed
seedm = int_m[idm,:];
#reverse interactions of seedmake
seedmrev = int_m[:,idm];

#Community composition
cid = [id;idm];

#Initialize the community matrix
c_m = vcat(seed,seedm);
crev_m = hcat(seedrev,seedmrev);

tmax = 100;

for t=1:tmax

  #Prior community structure
  c_mold = copy(c_m);
  crev_mold = copy(crev_m);

  #build a sublist of species that are trophically linked to anything in the community
  tlinksp = Array{Int64}(0);
  for i=1:length(cid)
    alink = find(x->x=='a',crev_m[:,i]);
    append!(tlinksp,alink);
  end
  #Randomize the list
  tlink = sample(tlinksp,length(tlinksp),replace=false)

  #Threshold check
  #run a while loop
  check = true;
  tic = 0;
  while check == true
    tic = tic + 1;
    #randomly select from the list
    draw = sample(tlinksp)

    dm = int_m[draw,:];
    dmrev = int_m[:,draw];

    #What things does this species need?
    draw_n = find(x->x=='n',dm);
    ln=length(draw_n);
    if ln>0
      #Are needs fulfilled to threshold?
      nperc = sum([in(draw_n[i],cid) for i=1:ln])/ln;
      if nperc >= n_thresh
        check = false
      end
    else
      check = false
    end
  end

  #Update the community
  push!(cid,draw);
  c_m = vcat(c_mold,int_m[draw,:]);
  crev_m = hcat(crev_mold,int_m[:,draw]);
