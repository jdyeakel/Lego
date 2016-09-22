function anime_sims_func(int_m,tp_m,tind_m,a_thresh,n_thresh,tmax)




  #Order of operations
  #1) Establish thresholds
  #2) Choose root species (must be primary producer)
  #3) Choose next species...
  #3a) determine and test 'a' and 'n' thresholds
  #3b) assess similarity to determine exclusion
  #3c) Pass/Fail
  #4) Cumulatively add species over time, where inclusion is determined by 1) thresholds, 2) exclusion


  # #Setting thresholds
  # n_thresh = 0.2;
  # a_thresh = 0.0;
  num_plays = length(int_m[1,:]);

  Slist = find(x->x=='n',diag(int_m));
  lS = length(Slist);

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
  #What things does the seed species need?
  idn = find(x->x=='n',seed);

  #Interactions of idseed
  seedm = int_m[idm,:];
  #reverse interactions of seedmake
  seedmrev = int_m[:,idm];

  #Community composition
  cid = [id;idm];

  #Initialize the community matrix
  new_play = 1+length(idm);
  c_m = Array{Char}(1+length(idm),num_play);
  crev_m = Array{Char}(num_play,1+length(idm));
  c_m[1,:] = seed;
  crev_m[:,1] = seedrev;
  for i=2:new_play
    c_m[i,:] = seedm;
    crev_m[:,i] = seedmrev;
  end

  #Store assembling community as a sparse matrix?
  com_sparse = Array{Char}(num_play,num_play);
  com_sparse[1:num_play^2]='0';
  for i=1:length(cid);
    com_sparse[cid[i],:] = c_m[i,:];
    com_sparse[:,cid[i]] = crev_m[:,i];
  end

  #Location on the species-only list?
  #The +1 accounts for the fact that the sun is included in species-only matrices (but not species-only list)
  spid = find(x->x==id,Slist)+1;
  #Update direct and indirect trophic interaction matrices
  com_tp = Array{Int64}(lS+1,lS+1)*0;
  com_tind = Array{Int64}(lS+1,lS+1)*0;
  #These matrices record species only
  com_tp[spid,:] = tp_m[spid,:];
  com_tp[:,spid] = tp_m[:,spid];
  com_tind[spid,:] = tind_m[spid,:];
  com_tind[:,spid] = tind_m[:,spid];


  c_t = Array{Char}(num_play,num_play,tmax);
  c_tp = Array{Int64}(lS+1,lS+1,tmax);
  c_tind = Array{Int64}(lS+1,lS+1,tmax);
  tout = Array{Int64}(1);

  time_tic = 0;
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


    did = Array{Char}(0);
    dseed = Array{Char}(0,num_play);
    dseedrev = Array{Char}(num_play,0);

    #Run this look IF tlink has elements
    #Threshold check
    #run a while loop to find the next 'colonizer'
    #If it fails, find another primary producer
    keepgoing = true;
    check = true;
    tic = 0;
    if length(tlink) > 0
      while check == true
        tic = tic + 1;

        ncheck = true;
        acheck = true;

        #randomly select from the list
        did = Array{Int64}(1);
        did[1] = tlink[tic];

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
      # print("Community is full at t=", t)
      break
    end

    #When nothing can colonize
    if keepgoing == false
      # print("Community is uninvadible at t=", t)
      break
    end

    #If we get here, the choice has 'passed' threshold analysis

    #Location on the species-only list?
    #The +1 accounts for the fact that the sun is included in species-only matrices (but not species-only list)
    spdid = find(x->x==did[1],Slist)+1;
    #Direct trophic interactions
    com_tp[spdid,:] = copy(tp_m[spdid,:]);
    com_tp[:,spdid] = copy(tp_m[:,spdid]);
    #Indirect trophic interactions
    com_tind[spdid,:] = copy(tind_m[spdid,:]);
    com_tind[:,spdid] = copy(tind_m[:,spdid]);

    #What does this species make?
    #Made things come along too!
    didm = find(x->x=='m',dseed);
    dm = int_m[didm,:];
    dmrev = int_m[:,didm];

    idnew = cat(1,did,didm);
    cmnew = vcat(dseed,dm);
    crevmnew = hcat(dseedrev,dmrev);
    #Update the sparse and trophic matrices
    for i=1:length(idnew);
      com_sparse[idnew[i],:] = copy(cmnew[i,:]);
      com_sparse[:,idnew[i]] = copy(crevmnew[:,i]);
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
    c_tp[:,:,t] = copy(com_tp);
    c_tind[:,:,t] = copy(com_tind);

    time_tic = time_tic + 1;

  end #end time loop

  return(time_tic, c_t, c_tp, c_tind, cid)

end
