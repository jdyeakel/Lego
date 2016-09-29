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
  check = true;
  while check == true
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
        #Needed things are in the community
        #Pass this test (false = pass)
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
      #Add sun as a 'consumable' within the community
      #This will allow true primary producers to colonize
      cid_wsun = cat(1,1,cid);
      aperc = sum([in(da[i],cid_wsun) for i=1:lda])/lda;
      # '>' because there MUST be at least one 'a' interaction
      if aperc > a_thresh
        #Assimilated things are in the community
        #Pass this test (false = pass)
        acheck = false;
      end
    else
      #This shouldn't happen based on constraints for tlink
      #But just for completeness - if you have no (a) links, then you do not pass this test
      acheck = true;
    end

    #If needs and assimilates are above thresholds, then stop while loop
    if ncheck == false && acheck == false
      check = false;
    end

  end





end
