function dynstructure(cid,cid_old)
    
    spcid = intersect(sp_v,cid);
    spcid_ind = indexin(spcid,[1;sp_v]);
    
    #Richness measures
    rich = length(cid);
    sprich = length(spcid);
    obrich = rich - sprich;
    
    #Turnover
    turnover = 1 - (length(intersect(cid,cid_old)) / length(union(cid,cid_old)));
    
    #Resource overlap
    res_overlap = roverlap(cid);
    #mean resource overlap
    mres_overlap = mean(res_overlap[!isnan.(res_overlap)]);
    
    return(rich,sprich,turnover,mres_overlap)
    
end