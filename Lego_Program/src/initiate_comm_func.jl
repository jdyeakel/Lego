function initiate_comm_func(int_m)

  #Species-only list
  Slist = find(x->x=='n',diag(int_m));
  lS = length(Slist);

  #Founding species
  #primary species
  prim_prod = find(x->x=='a',int_m[:,1]);
  id = rand(prim_prod);

  #update primary producer list
  deleteat!(prim_prod,find(x->x==id,prim_prod));
  #interactions of the seed species
  seed = copy(int_m[id,:]);
  #reverse interactions of the seed species
  seedrev = copy(int_m[:,id]);
  #What things does the seed species make?
  idm = find(x->x=='m',seed);

  #NOTE: Currently ignoring things that the initial colonizer 'needs'!
  #What things does the seed species need?
  idn = find(x->x=='n',seed);

  #Interactions of idseed
  seedm = int_m[idm,:];
  #reverse interactions of seedmake
  seedmrev = int_m[:,idm];

  #Community composition
  cid = [id;idm];

  #Initialize the community matrix
  #New players include the seed and its objects
  new_play = length(cid);
  c_m = Array{Char}(new_play,num_play);
  crev_m = Array{Char}(num_play,new_play);
  c_m[1,:] = seed;
  crev_m[:,1] = seedrev;
  for i=2:new_play
    c_m[i,:] = seedm[i-1,:];
    crev_m[:,i] = seedmrev[:,i-1];
  end

  # #Store assembling community as a sparse matrix?
  # com_sparse = Array{Char}(num_play,num_play);
  # com_sparse[1:num_play^2]='0';
  # for i=1:length(cid);
  #   com_sparse[cid[i],:] = c_m[i,:];
  #   com_sparse[:,cid[i]] = crev_m[:,i];
  # end

  #Location on the species-only list?
  #The +1 accounts for the fact that the sun is included in species-only matrices (but not species-only list)
  spid = find(x->x==id,Slist)+1;
  #Update direct and indirect trophic interaction matrices
  #Dimensions: number of species + the sun (1)
  com_tp = Array{Int64}(lS+1,lS+1)*0;
  com_tind = Array{Int64}(lS+1,lS+1)*0;
  #These matrices record species only
  com_tp[spid,:] = tp_m[spid,:];
  com_tp[:,spid] = tp_m[:,spid];
  com_tind[spid,:] = tind_m[spid,:];
  com_tind[:,spid] = tind_m[:,spid];

  return(cid, c_m, crev_m, com_tp, com_tind)

end
