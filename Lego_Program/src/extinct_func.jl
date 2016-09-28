function extinct_func(cid,c_m,crev_m,com_sparse,com_tp,com_tind)

  #Prior community structure
  cidold = copy(cid);
  c_mold = copy(c_m);
  crev_mold = copy(crev_m);

  num_com = length(cid);
  L = (sum(com_tp)/2);
  Slist = find(x->x=='n',diag(int_m));
  lS = length(Slist);

  #Determine the Predation and Competition Load for each species
  pred_vec = Array{Float64}(num_com);
  vuln_vec = Array{Float64}(num_com);
  comp_vec = Array{Float64}(num_com);
  for i=1:num_com
    intrev = int_m[:,cid[i]];
    inti = int_m[cid[i],:];
    #How many things ASSIMILATE it?
    pred_vec[i] = length(find(x->x=='a',intrev));
    vuln_vec[i] = (1/(L/lS))*pred_vec[i]
    #What is the max similarity?
    comp_load = Array{Float64}(num_com);
    for j=1:num_com
        intj = int_m[cid[j],:];
        seq1 = copy(inti);
        seq2 = copy(intj);
        comp_load[j] = sim_func(seq1,seq2);
    end
    #set itself to zero
    comp_load[i] = 0.0;
    comp_vec[i] = maximum(comp_load);
  end
  vulnscaled_vec = vuln_vec ./ maximum(vuln_vec);

  #Randomly select the species to eliminate
  esp = rand(cid);

  #Determine the cascade of consequent extinctions

  #1) What does it make? eliminate those things
  em = find(x->x=='m',int_m[esp,:]);


end
