function extinct_func(cid,c_m,crev_m,com_sparse,com_tp,com_tind)

  #Prior community structure
  cidold = copy(cid);
  c_mold = copy(c_m);
  crev_mold = copy(crev_m);

  #Number of links and species of prior community structure
  num_com = length(cid);
  L = (sum(com_tp)/2);
  Slist = find(x->x=='n',diag(int_m));
  lS = length(Slist);

  #Which cid-members are species?
  # spcid = Array{Int64}(0);
  # for i=1:num_com
  #   if in(cid[i],Slist) == true
  #     append!(spcid,cid[i]);
  #   end
  # end
  # num_comsp = length(spcid);

  #This is also encoded in the com_sparse matrix
  spcid = find(x->x=='n',diag(com_sparse));
  num_comsp = length(spcid);

  ########################
  # PRIMARY EXTINCTIONS
  ########################

  #Determine the Predation and Competition Load for each species
  pred_vec = zeros(num_comsp);
  vuln_vec = zeros(num_comsp);
  comp_vec = zeros(num_comsp);
  for i=1:num_comsp

    #Vector (species i interactions)
    inti = int_m[spcid[i],:];
    #Reverse vector (those interacting with species i)
    intrev = int_m[:,spcid[i]];

    #How many things ASSIMILATE it?
    #Found in the reverse vector
    pred_vec[i] = length(find(x->x=='a',intrev));
    vuln_vec[i] = (1/(L/lS))*pred_vec[i]

    #What is the max similarity?
    comp_load = Array{Float64}(num_comsp);
    for j=1:num_comsp
      intj = int_m[spcid[j],:];
      seq1 = copy(inti);
      seq2 = copy(intj);
      comp_load[j] = sim_func(seq1,seq2);
    end

    #set itself to zero
    comp_load[i] = 0.0;
    comp_vec[i] = maximum(comp_load);
  end
  #vulnscaled_vec = vuln_vec ./ maximum(vuln_vec);

  #Calculate the probability of extinction based on predation load and competitive load
  prext_pred = (1./(1.+exp(-0.2.*(pred_vec.-(L/lS)))));
  prext_comp = (1./(1.+exp(-10.*(comp_vec.-0.5))));
  prext = prext_pred.*prext_comp;

  #Determine extinction with a binomial draw based on the probability of extinction calculated above
  binext = Array{Int64}(num_comsp);
  for i = 1:num_comsp
    binext[i] = rand(Binomial(1,prext[i]));
  end

  #Select the species to eliminate
  esploc = find(x->x==1,binext);
  esp = spcid[esploc];
  lesp = length(esp);
  #Eliminate those species from the spcid list
  #Note: esploc must be sorted for the deleteat! function
  #deleteat!(spcid,sort(esploc));

  #Record the species-only eliminations
  #i.e. this EXCLUDES made things for updating trophic nets
  esponly = copy(esp);

  #1) What does it make? eliminate those things IF nothing else present makes them either
  for i=1:lesp
    #The made objects of the species being deleted
    madeobjects = find(x->x=='m',int_m[esp[i],:]);
    #Only do this loop if things are made
    if length(madeobjects) > 0
      #Iterate across each thing that is MADE
      #what things make this?
      for k=1:length(madeobjects)
        made = madeobjects[k];
        makers = find(x->x=='n',int_m[made,:]);
        #Is there anything that makes this in the community that is NOT going extinct during this timestep?
        checkmade = Array{Bool}(length(makers));
        for j=1:length(makers)
          if in(makers[j],cid) && in(makers[j],esponly)==false
            #There is something else that makes it
            checkmade[j] = true;
          else
            #There is NOT something else that makes it
            checkmade[j] = false;
          end
        end
        #If there are NO makers in the community that are remaining
        #then add the made thing to the eliminate list
        if all(x->x==false,checkmade)
          append!(esp,made);
          append!(esploc,find(x->x==made,cid));
        end
      end
    end
  end

  #UPDATING
  #Update CID by eliminating esp
  for i=1:lesp
    del_loc = find(x->x==esp[i],cid);
    deleteat!(cid,del_loc);
  end

  #Update c_m, crev_m,
  #Species interactions
  c_m = copy(int_m[cid,:]);
  #Interactions ON each species
  crev_m = copy(int_m[:,cid]);

  #Update interaction matrices from the esponly vector
  #This is the vector of species-only deletions

  #NOTE: This doesn't work bc species id doesn't correlate with the row/columns of the trophic matrices
  for i=1:length(esponly)
    com_tp[esponly[i],:] = zeros(Int64,(lS+1));
    com_tp[:,esponly[i]] = zeros(Int64,(lS+1));
    com_tind[esponly[i],:] = zeros(Int64,(lS+1));
    com_tind[:,esponly[i]] = zeros(Int64,(lS+1));
  end

  #Update com_sparse?
  com_sparse[esp,:] = '0';
  com_sparse[:,esp] = '0';
  #Replace diagonal of com_sparse to original from int_m
  #NOTE: THIS PART DOESNT WORK YET
  diag(com_sparse) = copy(diag(int_m));
  #NOTE: to keep it consistent, change 'i' to '0' though I don't know why we have to do that
  diag(com_sparse)[find(x->x=='i')]='0';

  #Update the number of species in the community
  num_com = length(cid);
  num_comsp = length(spcid);

  ########################
  # SECONDARY EXTINCTIONS
  ########################

  #Check ASSIMILATION AND NEED THRESHOLDS after esp has been eliminated
  check = true;
  while check == true
    ncheck = trues(num_comsp);
    acheck = trues(num_comsp);
    ancheck = trues(num_comsp);
    #Cycle through the SPECIES in the community (not made things)
    for i=1:num_comsp
      #What species?
      did = spcid[i];
      dseed = int_m[did,:];
      #CHECKING NEEDS
      #What things does this species need?
      dseedneed = copy(dseed);
      dseedneed[did] = '0';
      dn = find(x->x=='n',dseedneed);
      ldn=length(dn);
      nperc = Array{Float64}(1);
      if ldn>0
        #Are needs fulfilled to threshold?
        needs = Array{Bool}(ldn);
        for j=1:ldn
          needs[j] = in(dn[j],cid);
        end
        nperc = sum(needs)/ldn;
        if nperc >= n_thresh
          #Needed things are in the community
          #Pass this test (false = pass)
          ncheck[i] = false;
        end
      else
        #Chosen immigrant needs nothing
        #Move on to next step (false = pass)
        ncheck[i] = false;
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
        prey = Array{Bool}(lda);
        for j=1:lda
          prey[j] = in(da[j],cid_wsun);
        end
        aperc=sum(prey)/lda;
        # '>' because there MUST be at least one 'a' interaction
        if aperc > a_thresh
          #Assimilated things are in the community
          #Pass this test (false = pass)
          acheck[i] = false;
        end
      else
        #This shouldn't happen based on constraints for tlink
        #But just for completeness - if you have no (a) links, then you do not pass this test
        acheck[i] = true;
      end

      #If needs and assimilates are above thresholds, then stop while loop
      if ncheck[i] == false && acheck[i] == false
        #if you have pass your needs and assimilates threshold, you live
        #false = pass
        ancheck[i] = false;
      end

    end #End for loop

    if all(i->i==false,ancheck)
      check = false; #this will break the while loop
    end

    #Find & impliment secondary xtinctions
    find(x->x==true,ancheck)



  end #End while loop





end #End function
