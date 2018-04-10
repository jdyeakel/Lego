function roverlap(cid,sp_v,a_b)
    
    # cids = sort(cid);
    
    spcid = intersect(sp_v,cid);
    l_sp = length(spcid);
    
    
    
    #Species consuming each resource (species and objects)
    preds = vec(sum(a_b[cid,cid],1));
    res = vec(sum(a_b[spcid,cid],2));
    
    
    # resource overlap only tallied for species, but resources include objects
    
    # res_overlap_old = Array{Float64}(l_sp);
    # for i=1:length(spcid)
    #     prey = find(isodd,a_b[spcid[i],cid]);
    #     res_overlap_old[i] = mean((preds[prey]-1)./(l_sp-1));
    # end
    
    #Alternative calculate
    #10X faster
    res_overlap = (((a_b[spcid,cid]*preds).-res)/(l_sp-1))./res;
    
    #This is the same thing, but the above code is 10x faster!
    # res_overlap = (((a_b[spcid,cid]*preds).-res)./res)/(l_sp-1);
    
    #Obligate primary producers will be NaN
    #They should be discluded from the measurement so we will keep them NaN
    
    return(res_overlap)
    
end