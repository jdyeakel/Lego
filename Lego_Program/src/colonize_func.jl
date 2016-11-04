function colonize_func(int_m,tp_m,tind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind)

  num_play = length(int_m[1,:]);

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

  #build a sublist of species that are trophically linked to anything in the community - this will involve only species
  #Start with the primary producers
  tlink_full = copy(prim_prod);
  for i=1:length(cidold)
    alink = find(x->x=='a',crev_m[:,i]);
    append!(tlink_full,alink);
  end

  #Get rid of duplicates
  tlink_unique = unique(tlink_full);

  #Eliminate anything that is already there
  #1) Build a to-delete list
  lt = length(tlink_unique);
  torm = Array{Int64}(0);
  for i=1:lt
    if in(tlink_unique[i],cidold)
      #The positions of things to delete
      push!(torm,i);
    end
  end
  #2) delete those already there
  deleteat!(tlink_unique,torm);


  #Randomize the list
  #We will walk through the list to find a potential colonizer, however we don't want the order to bias selection.
  tlink = sample(tlink_unique,length(tlink_unique),replace=false);


  # did = Array{Char}(0);
  # dseed = Array{Char}(0,num_play);
  # dseedrev = Array{Char}(num_play,0);

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
      did = copy(tlink[tictoc]);

      dseed = copy(int_m[did,:]);
      dseedrev = copy(int_m[:,did]);

      #CHECKING NEEDS
      #What things does this species need?
      dseedneed = copy(dseed);
      dseedneed[did] = '0';
      dn = find(x->x=='n',dseedneed);
      ldn=length(dn);
      nperc = zeros(Float64,1);
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
      aperc = zeros(Float64,1);
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

      #If needs and assimilates are above thresholds, then the colonizer passes and we can stop while loop
      #False = pass
      if ncheck == false && acheck == false
        check = false;
      end

      #Nothing can colonize?
      #I.e. check is true (so there was a fail in the need or assimilate test) and we have reached the end of the list
      if check == true && tictoc == length(tlink)
        #To stop the while loop
        check = false;
        #To stop the colonization module
        keepgoing = false;
      end

    end #End while loop

  else
    #If nothing in the template eats anything in the community, then the community is full and we stop the module
    status = "full";
    #println("Community is trophically disconnected at t=", t)
    return(status,cid,c_m,crev_m,com_tp,com_tind);
  end

  #If we have reached the end of the list of potential colonizers and nothing passes the need or assimilate tests, the community is full and we stop the module
  if keepgoing == false
    status = "full";
    #println("Community is uninvadible at t=", t)
    return(status,cid,c_m,crev_m,com_tp,com_tind);
  else
    #If we get here, the choice has 'passed' threshold analysis
    #The community is still open, and we must now update the community to reflect the new added colonizers and the things that they make

    status = "open";

    #NOTE: I think there is an error in updating the trophic (direct and indirect) matrices

    #Updating the community trophic matrices
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

    # #Update the sparse matrices
    # for i=1:lnew;
    #   com_sparse[idnew[i],:] = copy(cmnew[i,:]);
    #   com_sparse[:,idnew[i]] = copy(crevmnew[:,i]);
    # end

    #Update the community
    cid_update = cat(1,cidold,idnew);
    c_m_update = vcat(c_mold,cmnew);
    crev_m_update = hcat(crev_mold,crevmnew);

    #Update the updated variables to in-play variables
    cid = copy(cid_update);
    c_m = copy(c_m_update);
    crev_m = copy(crev_m_update);

    return(status,cid,c_m,crev_m,com_tp,com_tind);
  end

end
