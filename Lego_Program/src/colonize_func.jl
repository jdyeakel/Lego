function colonize_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind)

  num_play = length(int_m[1,:]);

  status = "open";

  #Species-only list
  Slist = find(x->x=='n',diag(int_m));
  lS = length(Slist);

  #Prior community structure
  cidold = copy(cid);
  c_mold = copy(c_m);
  crev_mold = copy(crev_m);
  
  # # Make a list of species-only in community
  # intcom = int_m[cid,cid];
  # spcid = find(x->x=='n',diag(intcom));
  sponly = copy(cid);
  dl = Array{Int64}(0);
  for i=1:length(cid)
    if int_m[sponly[i],sponly[i]] == 'i'
      push!(dl,i);
    end
  end
  deleteat!(sponly,dl);
  lsp = length(sponly);
  spcid = copy(sponly);
  
  
  #List of primary producers
  prim_prod = find(x->x=='a',int_m[:,1]);

  #build a sublist of species that are trophically linked to anything in the community - this will involve only species
  #Start with the primary producers
  tlink_full = copy(prim_prod);
  for i=1:length(cidold)
    #What EATS cid[i]? Where crev_m[:,i] is 'a'
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


  did = Array{Char}(0);
  dseed = Array{Char}(0,num_play);
  dseedrev = Array{Char}(num_play,0);

  #Run this look IF tlink has elements
  #Threshold check
  #run a while loop to find the next 'colonizer'
  #If it fails, find another primary producer
  
  #This is an exhaustive search because we want to count the NUMBER of potential colonizers
  #In contrast, the colonizesingle_func.jl does not perform an exhaustive search, and only searches until a colonizer is found, so is faster.
  #potcol records the potential colonizers (does not count objects)
  potcol = Array{Int64}(0);
  
  if length(tlink) > 0
    for i=1:length(tlink)
    #while check == true
      #tictoc = tictoc + 1;

      ncheck = true;
      acheck = true;

      #randomly select from the list
      did = copy(tlink[i]);

      dseed = copy(int_m[did,:]);
      dseedrev = copy(int_m[:,did]);

      #CHECKING NEEDS
      #What things does this species need?
      dseedneed = copy(dseed);
      #the diagonal 'need' just identifies it as a species (not an object), and does not indicate an interaction need. Set this to '0' to discount it.
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
        #Are assimilates fulfilled to threshold?
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
        #check = false;
        #Colonizer _did CAN fit into the community, so add it to potcol list
        append!(potcol,did)
      end

      #Nothing can colonize?
      #I.e. check is true (so there was a fail in the need or assimilate test) and we have reached the end of the list
      # if check == true && tictoc == length(tlink)
      #   #To stop the while loop
      #   check = false;
      #   #To stop the colonization module
      #   keepgoing = false;
      # end

    end #End for loop

  else
    #If nothing in the template eats anything in the community, then the community is full and we stop the module
    status = "full";
    #println("Community is trophically disconnected at t=", t)
    return(status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,potcol);
  end
  
  
  #If there are no colonizers that PASS, end the colonization function
  if length(potcol) == 0
    status = "full";
    #println("Community is uninvadible at t=", t)
    return(status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,potcol);
  else
    #If we get here, the choice has 'passed' threshold analysis
    #The community is still open, and we must now update the community to reflect the new added colonizers and the things that they make

    status = "open";
    
    #Choose the colonizer from the list of potcol
    #Restablish its list of interactions
    did = rand(potcol);
    dseed = copy(int_m[did,:]);
    dseedrev = copy(int_m[:,did]);
    
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
      #Delete the focal species that is to be added to the community
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


    #Update trophic and mutualistic interactions

    #1) make a list of species-only in community
    sponly = copy(cid);
    dl = Array{Int64}(0);
    for i=1:length(cid)
      if int_m[sponly[i],sponly[i]] == 'i'
        push!(dl,i);
      end
    end
    deleteat!(sponly,dl);
    lsp = length(sponly);
    spcid = copy(sponly);

    #2) Find the location of species on trophic matrix
    #This will be used to locate the correct species on trophic matrices
    t_loc = zeros(Int64,lsp);
    for i=1:lsp
      sploc = find(x->x==spcid[i],Slist);
      if length(sploc)>1
        print("Duplication Error")
      end
      #+1 because the 1st row/column of the trophic matrices is the sun
      t_loc[i] = sploc[1] + 1;
    end

    #Attach the primary producer
    t_loc = copy([1;t_loc]);

    #3) Sync the community trophic matrix with the template trophic matrix FOR ONLY SPECIES/INTERACTIONS THAT ARE PRESENT
    #Direct trophic interactions
    com_tp[t_loc,t_loc] = copy(tp_m[t_loc,t_loc]);
    #Indirect trophic interactions
    com_tind[t_loc,t_loc] = copy(tind_m[t_loc,t_loc]);
    #Mutualistic interactions
    com_mp[t_loc,t_loc] = copy(mp_m[t_loc,t_loc]);
    #Indirect mutualistic interactions
    com_mind[t_loc,t_loc] = copy(mind_m[t_loc,t_loc]);

    return(status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,potcol);
  end

end
