function potcol(sp_v,int_id,cid,a_b,n_b0,athresh,nthresh)
    #Which are species?
    spcid = intersect(sp_v,cid);
    #Which are objects?
    ocid = setdiff(cid,spcid);
    
    
    #COUNT POTENTIAL COLONIZERS
    trophiclinked = setdiff(int_id[(sum(a_b[:,[1;cid]],2) .> 0)[:,1]],cid);
    #For each trophiclinked, count number of assimilate and need interactions in system
    #Determine in the proportion that already exists is >= the threshold
    a_fill = ((sum(a_b[trophiclinked,[1;cid]],2)./sum(a_b[trophiclinked,:],2)) .>= athresh)[:,1];
    prop_n = sum(n_b0[trophiclinked,[1;cid]],2)./sum(n_b0[trophiclinked,:],2);
    #If there are no 'need' interactions, proportion filled is always 1, and will always pass
    prop_n[isnan.(prop_n)] = 1;
    n_fill = (prop_n .>= nthresh)[:,1];
    #Build a list of those that pass each individual test
    a_pass = trophiclinked[a_fill];
    n_pass = trophiclinked[n_fill];
    #Build a list of those that pass both tests (potential colonizers)
    col = intersect(a_pass,n_pass);
    #Count the number that pass
    lcol = length(col);
    
    return col
end
