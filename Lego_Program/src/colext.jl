function colext(int_m,cid,a_thresh,n_thresh)
    
    #############
    # Preamble... could take this out of the function
    
    #Boolian matrices
    a_b = (int_m .== 'a')*1; 
    n_b = (int_m .== 'n')*1;
    i_b = (int_m .== 'i')*1;
    m_b = (int_m .== 'm')*1;
    
    #copy of need binary matrix with diag = 0
    n_b0 = copy(n_b);
    n_b0[diagind(n_b0)]=0;

    #Vector and length of species IDs
    sp_v = find(isodd,diag(n_b));
    l_sp = length(sp_v);
    
    int_id = collect(1:l_intm);

    
    ######################
    #CURRENT COMMUNITY
    ######################
    
    #ID's of agents in the community
    cid_old = copy(cid);
    
    #IDs of species only
    spcid = intersect(sp_v,cid_old);
    
    ######################
    #COLONIZATION MODULE
    ######################
    
    #which species have an 'a'-link to anyone in the community? (this should include primary producers and those that are trophically linded to objects: [1 cid])
    trophiclinked = setdiff(int_id[(sum(a_b[:,[1;cid]],2) .> 0)[:,1]],cid);
    
    #Define the status of the community as open to colonizers
    status = "open";
    
    #For each trophiclinked, count number of assimilate and need interactions in system
    #Determine in the proportion that already exists is >= the threshold
    a_fill = ((sum(a_b[trophiclinked,[1;cid]],2)./sum(a_b[trophiclinked,:],2)) .>= a_thresh)[:,1];
    prop_n = sum(n_b0[trophiclinked,[1;cid]],2)./sum(n_b0[trophiclinked,:],2);
    #If there are no 'need' interactions, proportion filled is always 1, and will always pass
    prop_n[isnan(prop_n)] = 1;
    n_fill = (prop_n .>= n_thresh)[:,1];
    
    #Build a list of those that pass each individual test
    a_pass = trophiclinked[a_fill];
    n_pass = trophiclinked[n_fill];
    
    #Build a list of those that pass both tests (potential colonizers)
    pot_col = intersect(a_pass,n_pass);
    
    #Count the number that pass
    lpot_col = length(pot_col);
    
    #Choose a colonizer at random
    c_sp = rand(pot_col);
    
    #select made objects as well
    c_ob_full = find(isodd,m_b[c_sp,:]);
    
    #Only bring in objects that aren't already in the community
    c_ob = setdiff(c_ob_full,cid_old);
    
    colonize = [c_sp;c_ob];
    
    if (length(colonize) == 0)
        status = "closed";
    end
    
    #Update community
    cid = [cid_old;colonize];
    
    ######################
    #EXTINCTION MODULE
    ######################
    
    #run the extinction module or colonization by itself?
    if extinctions == true
        
        cid_old = copy(cid);
        
        #which are species?
        spcid = intersect(sp_v,cid_old);
        #which are objects?
        ocid = setdiff(cid_old,spcid);
        
        #1) determine predation load for each species
        num_preds = sum(a_b[spcid,spcid],2)[:,1];
        
        #2) Caluculate extinction probability
        avgk = 8; # avgk = gL/gS;
        prext_pred = (1./(1+exp(-0.5.*(num_preds.-(trophicload*avgk)))));
        
        #3) Draw extinctions
        binext = rand.(Binomial.(1,prext_pred));
        
        num_ext1 = sum(binext);
        
        #Only go through this round if there are any extinctions
        if num_ext1 > 0
            
            #####################
            # Primary extinctions
            #####################
            
            #Which species are going extinct?
            e_sp = spcid[find(isodd,binext)];
            le_sp = length(e_sp);
            
            #Delete species and their UNIQUELY MADE objects
            #what objects are uniquely made by esp?
            # e_ob = ocid[(sum(m_b[e_sp,ocid],1).==1)[1,:]];
            e_ob = ocid[find(isodd,sum(m_b[e_sp,ocid],1))]
            le_ob = length(e_ob);
            
            prim_extinct = [e_sp;e_ob];
            
            #Update community
            cid = setdiff(cid_old,prim_extinct);
            
        
            
            #######################
            # Secondary extinctions
            #######################
            
            #Only do this if there are still members in the community
            if length(cid) > 0
            
                num_ext2 = 0;
                
                #The assimilate/need test has to pass (an_test == true) to avoid any more secondary extinctions
                an_test = false;
                
                while an_test == false
                    #Which are species?
                    spcid = intersect(sp_v,cid);
                    #Which are objects?
                    ocid = setdiff(cid,spcid);
                
                    #Re-calculate a_thresh and n_thresh
                    a_fill = ((sum(a_b[spcid,[1;cid]],2)./sum(a_b[spcid,:],2)) .>= a_thresh)[:,1];
                    prop_n = sum(n_b0[spcid,[1;cid]],2)./sum(n_b0[spcid,:],2);
                    #If there are no 'need' interactions, proportion filled is always 1, and will always pass
                    prop_n[isnan(prop_n)] = 1;
                    n_fill = (prop_n .>= n_thresh)[:,1];
                    
                    #Build a list of those that pass each individual test
                    a_pass = spcid[a_fill];
                    n_pass = spcid[n_fill];
                    
                    #which species pass?
                    survivors = intersect(a_pass,n_pass);
                    
                    #If there are secondary extinctions...
                    if length(survivors) != length(spcid)
                    
                        #Species that secondarily don't pass the a-n test
                        e_sp2 = setdiff(spcid,survivors);
                        
                        #Number of secondary extinctions
                        le_sp2 = length(e_sp2);
                        
                        num_ext2 += le_sp2;
                        
                        #Delete species and their UNIQUELY MADE objects
                        #what objects are uniquely made by esp?
                        e_ob2 = ocid[(sum(m_b[e_sp2,ocid],1).==1)[1,:]];
                        le_ob2 = length(e_ob);
                        
                        sec_extinct = [e_sp;e_ob];
                        
                        #Update community
                        cid = setdiff(cid,sec_extinct);
                    
                    else
                        
                        #If there are no survivors, test is passed; exit loop
                        an_test = true;
                    
                    end
                    
                end #end secondary extinction module
            
            end #end if length(cid) > 0
        
        end #end extinctions module
        
    end #end (if extinctions==true) module
    
    return(
    cid,
    lpot_col,
    status,
    num_ext1,
    num_ext2
    )
    
end #end colext function
