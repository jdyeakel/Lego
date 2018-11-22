function roverlap(cid,sp_v,a_b,n_b0)
    
    # cids = sort(cid);
    
    spcid = intersect(sp_v,cid);
    l_sp = length(spcid);
    
    
    #TROPHIC OVERLAP
    #Species consuming each resource (species and objects)
    # preds = vec(sum(a_b[cid,cid],1));
    # res = vec(sum(a_b[spcid,cid],2));
    # res_overlap = (((a_b[spcid,cid]*preds).-res)/(l_sp-1))./res;
    
    
    
    #OR IF WE WANT TO CALCULATE TOTAL RESOURCE (EAT + NEED) OVERLAP
    #Species consuming or needing each resource (species and objects)
    A = a_b[spcid,cid] .+ n_b0[spcid,cid];
    
    # users = vec(sum(a_b[cid,cid],1)) .+ vec(sum(n_b0[cid,cid],1));
    users = vec(sum(a_b[spcid,cid],dims=1)) .+ vec(sum(n_b0[spcid,cid],dims=1));
    
    used = vec(sum(a_b[spcid,cid],dims=2)) .+ vec(sum(n_b0[spcid,cid],dims=2));
    #Proporitonal overlap of resources between species
    res_overlap = (((((a_b[spcid,cid] .+ n_b0[spcid,cid]))*users) .- used)/(length(spcid)-1)) ./ used;
    
    user_overlap = (((((a_b[spcid,cid] .+ n_b0[spcid,cid]))'*used) .- users)/(length(spcid)-1)) ./ users;
    #Save only overlap in objects
    user_overlap = user_overlap[1:l_sp];
    
    #species with no users have value nan; set this equal to 0
    user_overlap[findall(isnan,user_overlap)] .= 0;
    
    #Obligate primary producers will be NaN
    #They should be discluded from the measurement so we will keep them NaN
    
    
    return(res_overlap,user_overlap)
    
end
