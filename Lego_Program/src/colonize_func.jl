function colonize_func(a_thresh,n_thresh,cid,c_m,crev_m,com_sparse,com_tp,com_tind)

  status = "open";

  #Species-only list
  Slist = find(x->x=='n',diag(int_m));
  lS = length(Slist);

  #Prior community structure
  cidold = copy(cid);
  c_mold = copy(c_m);
  crev_mold = copy(crev_m);

  #List of primary producers
  prim_prod = find(x->x=='a',int_m[:,1]);

  #build a sublist of species that are trophically linked to anything in the community
  #Start with the primary producers
  tlink_full = copy(prim_prod);
  for i=1:length(cidold)
    alink = find(x->x=='a',crev_m[:,i]);
    append!(tlink_full,alink);
  end

  #Get rid of duplicate
  tlink_unique = unique(tlink_full);

  #Eliminate anything that is already there
  lt = length(tlink_unique);
  torm = Array{Int64}(0);
  for i=1:lt
    if in(tlink_unique[i],cidold)
      push!(torm,i);
    end
  end
  deleteat!(tlink_unique,torm);


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
  tictoc = 0;
  if length(tlink) > 0
    while check == true
      tictoc = tictoc + 1;

      ncheck = true;
      acheck = true;

      #randomly select from the list
      did = tlink[tictoc];

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

      #Nothing can colonize?
      if check == true && tictoc == length(tlink)
        #to stop the while loop
        check = false;
        keepgoing = false;
      end

    end #End while loop

  else
    status = "full";
    #println("Community is trophically disconnected at t=", t)
    return(status,cid,c_m,crev_m,com_sparse,com_tp,com_tind);
  end

  #When nothing can colonize
  if keepgoing == false
    status = "full";
    #println("Community is uninvadible at t=", t)
    return(status,cid,c_m,crev_m,com_sparse,com_tp,com_tind);
  else

    #If we get here, the choice has 'passed' threshold analysis
    status = "open";


    #Location on the species-only list?
    #The +1 accounts for the fact that the sun is included in species-only matrices (but not species-only list)
    spdid = find(x->x==did,Slist)+1;
    #Direct trophic interactions
    com_tp[spdid,:] = copy(tp_m[spdid,:]);
    com_tp[:,spdid] = copy(tp_m[:,spdid]);
    #Indirect trophic interactions
    com_tind[spdid,:] = copy(tind_m[spdid,:]);
    com_tind[:,spdid] = copy(tind_m[:,spdid]);

    #What does this species make?
    #Made things come along too!
    didm = find(x->x=='m',dseed);

    ldm = length(didm);

    #Only add on made things that aren't already made by something else in the community
    maderm = Array{Int64}(0);
    for i=1:ldm
      made = copy(didm[i]);
      #What other things make this?
      makers = find(x->x=='m',int_m[:,made]);
      #Delete the one that is to be added to the community
      deleteat!(makers,find(x->x==did,makers));
      #Are any of the other makers in the community?
      #If so, delete the made object, as its being added is a duplication
      for j=1:length(makers)
        if in(makers[j],cid)
          append!(maderm,i);
        end
      end
    end
    maderm_unique = unique(maderm);
    deleteat!(didm,sort(maderm_unique));

    #Recalculate ldm
    ldm = length(didm);

    dm = int_m[didm,:];
    dmrev = int_m[:,didm];

    idnew = cat(1,did,didm);
    lnew = length(idnew);
    #ldm + 1 to account for made things + maker
    cmnew = Array{Char}(lnew,num_play);
    crevmnew = Array{Char}(num_play,lnew);
    cmnew[1,:] = dseed;
    crevmnew[:,1] = dseedrev;
    for i=2:lnew
      cmnew[i,:] = dm[i-1,:];
      crevmnew[:,i] = dmrev[:,i-1];
    end
    # cmnew = vcat(dseed,dm);
    # crevmnew = hcat(dseedrev,dmrev);

    #Update the sparse and trophic matrices
    for i=1:lnew;
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

    return(status,cid,c_m,crev_m,com_sparse,com_tp,com_tind);
  end

end
