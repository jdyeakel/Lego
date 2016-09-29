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
  #vulnscaled_vec = vuln_vec ./ maximum(vuln_vec);

  #Calculate the probability of extinction
  prext_pred = (1./(1.+exp(-0.2.*(pred_vec.-(L/lS)))));
  prext_comp = (1./(1.+exp(-10.*(comp_vec.-0.5))));
  prext = prext_pred.*prext_comp;

  #Determine extinction with a binomial draw
  binext = Array{Int64}(num_com);
  for i = 1:num_com
    binext[i] = rand(Binomial(1,prext[i]));
  end

  #Randomly select the species to eliminate
  esploc = find(x->x==1,binext);
  esp = cid[esploc];

  #Determine the cascade of consequent extinctions

  #1) What does it make? eliminate those things
  for i=1:length(esp)
    made = find(x->x=='m',int_m[esp[i],:]);
    append!(esp,made);
    for j=1:length(made)
      append!(esploc,find(x->x==made[j],cid));
    end
  end

  #UPDATING
  #Update CID list by eliminating esp
  deleteat!(cid,sort(esploc));

  #Update c_m, crev_m,
  c_m = copy(int_m[cid,:]);
  crev_m = copy(int_m[:,cid]);

  #Update interaction matrices
  # com_sparse[esp,:] = Array{Int64}(num_play)*0;
  # com_sparse[:,esp] = Array{Int64}(num_play)*0;
  for i=1:length(esp)
    com_tp[esp[i],:] = Array{Int64}(lS+1)*0;
    com_tp[:,esp[i]] = Array{Int64}(lS+1)*0;
    com_tind[esp[i],:] = Array{Int64}(lS+1)*0;
    com_tind[:,esp[i]] = Array{Int64}(lS+1)*0;
  end

  #Check ASSIMILATION AND NEED THRESHOLDS after esp has been eliminated
  

end
