function roverlap(cid)
    
    # cids = sort(cid);
    
    spcid = intersect(sp_v,cid);
    l_sp = length(spcid);
    
    res_overlap = Array{Float64}(l_sp);
    
    #Species consuming each resource (species and objects)
    preds = vec(sum(a_b[cid,cid],1));
    
    #resource overlap only tallied for species, but resources include objects
    for i=1:length(spcid)
        prey = find(isodd,a_b[spcid[i],cid]);
        res_overlap[i] = mean((preds[prey]-1)./(l_sp-1));
    end
    #Obligate primary producers will be NaN
    #They should be discluded from the measurement so we will keep them NaN
    
    return(res_overlap)
    
end