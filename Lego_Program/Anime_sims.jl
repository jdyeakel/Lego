using Distributions
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")

num_play = 20

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
#update primary producer list
deleteat!(prim_prod,find(x->x==id,prim_prod));
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

#Store assembling community as a sparse matrix?
com_sparse = Array{Char}(num_play,num_play);
com_sparse[1:num_play^2]='0';
for i=1:length(cid);
  com_sparse[cid[i],:] = c_m[i,:];
  com_sparse[:,cid[i]] = crev_m[:,i];
end
c_t = Array{Char}(num_play,num_play,tmax);

for t=1:tmax

  #Prior community structure
  cidold = copy(cid);
  c_mold = copy(c_m);
  crev_mold = copy(crev_m);

  #build a sublist of species that are trophically linked to anything in the community
  #Start with the primary producers
  tlink_full = copy(prim_prod);
  for i=1:length(cidold)
    alink = find(x->x=='a',crev_m[:,i]);
    append!(tlink_full,alink);
  end
  #Eliminate anything that is already there
  lt = length(tlink_full);
  torm = Array{Int64}(0);
  for i=1:lt
    if in(tlink_full[i],cidold)
      push!(torm,i);
    end
  end
  deleteat!(tlink_full,torm);
  tlink_unique = unique(tlink_full);

  #Randomize the list
  tlink = sample(tlink_unique,length(tlink_unique),replace=false);

  #Run this look IF tlink has elements
  if length(tlink) > 0
    #Threshold check
    #run a while loop to find the next 'colonizer'
    #If it fails, find another primary producer
    keepgoing = true;
    check = true;
    tic = 0;
    while check == true
      tic = tic + 1;

      ncheck = true;
      acheck = true;

      #randomly select from the list
      did = tlink[tic];

      dseed = int_m[did,:];
      dseedrev = int_m[:,did];

      #CHECKING NEEDS
      #What things does this species need?
      dseedneed = copy(dseed);
      dseedneed[did] = '0';
      dn = find(x->x=='n',dseedneed);
      ldn=length(dn);
      nperc = Array{Float64}(1);
      if ldn>0
        #Are needs fulfilled to threshold?
        nperc = sum([in(dn[i],cid) for i=1:ldn])/ldn;
        if nperc >= n_thresh
          ncheck = false;
        end
      else
        #Chosen immigrant needs nothing
        #Move on to next step (false = pass)
        ncheck = false;
      end

      #CHECKING ASSIMILATES
      #What things does this species Assimilate?
      dseedass = copy(dseed);
      da = find(x->x=='a',dseedass);
      lda=length(da);
      aperc = Array{Float64}(1);
      if lda>0
        #Are needs fulfilled to threshold?
        aperc = sum([in(da[i],cid) for i=1:lda])/lda;
        if aperc >= a_thresh
          acheck = false;
        end
      else
        #This shouldn't happen based on constraints for tlink
        #But just for completeness
        acheck = false;
      end

      #If needs and assimilates are above thresholds, then stop while loop
      if ncheck == false && acheck == false
        check = false;
      end

      #Nothing can colonize?
      if check == true && tic == length(tlink)
        #to stop the while loop
        check = false;
        keepgoing = false;
      end

    end

  else
    print("Community is full at t=", t)
    break
  end

  #When nothing can colonize
  if keepgoing == false
    print("Community is uninvadible at t=", t)
    break
  end

  #If we get here, the choice has 'passed' threshold analysis

  #What does this species make?
  #Made things come along too!
  didm = find(x->x=='m',dseed);
  dm = int_m[didm,:];
  dmrev = int_m[:,didm];

  idnew = cat(1,did,didm);
  cmnew = vcat(dseed,dm);
  crevmnew = hcat(dseedrev,dmrev);
  #Update the sparse matrix
  for i=1:length(idnew);
    com_sparse[idnew[i],:] = cmnew[i,:];
    com_sparse[:,idnew[i]] = crevmnew[:,i];
  end

  #Update the community
  cid_update = cat(1,cidold,idnew);
  c_m_update = vcat(c_mold,cmnew);
  crev_m_update = hcat(crev_mold,crevmnew);

  #Update the updated variables to in-play variables
  cid = copy(cid_update);
  c_m = copy(c_m_update);
  crev_m = copy(crev_m_update);

  #Save image of the community through time
  c_t[:,:,t] = copy(com_sparse);

end #end time loop
