function initiate_comm_func(int_m,tp_m,tind_m)

  num_play = length(int_m[1,:]);

  #Species-only list
  Slist = find(x->x=='n',diag(int_m));
  lS = length(Slist);

  #Founding species
  #Find which species are primary producers
  prim_prod = find(x->x=='a',int_m[:,1]);
  lprim = length(prim_prod);
  #How many things do the primary producers need?
  pp_need = Array{Bool}(lprim);
  for i=1:lprim
    pp_need[i] = length(find(x->x=='n',int_m[prim_prod[i],:])) <2;
  end

  #Which primary producers are indendent?
  #Note: Currently species 2 is built as being independent, so there should always be at least one capable species
  prim_prod_indep = prim_prod[find(x->x==true,pp_need)];

  #Randomly choose a primary producer to seed the community
  id = rand(prim_prod_indep);

  # #update primary producer list
  # deleteat!(prim_prod_indep,find(x->x==id,prim_prod_indep));

  #interactions of the seed species
  seed = copy(int_m[id,:]);

  #reverse interactions of the seed species
  seedrev = copy(int_m[:,id]);

  #What things does the seed species make?
  idm = find(x->x=='m',seed);

  #Interactions of things made by seed
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

  #BUG: There is an error in updating the trophic (direct and indirect) matrices

  #Update direct and indirect trophic interaction matrices
  #Dimensions: number of species + the sun (1)
  com_tp = zeros(Int64,lS+1,lS+1);
  com_tind = zeros(Int64,lS+1,lS+1);
  #These matrices record species only
  com_tp[spid,:] = tp_m[spid,:];
  com_tp[:,spid] = tp_m[:,spid];
  com_tind[spid,:] = tind_m[spid,:];
  com_tind[:,spid] = tind_m[:,spid];

  return(cid, c_m, crev_m, com_tp, com_tind)

end
