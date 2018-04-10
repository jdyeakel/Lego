function dynstructure(cid,cid_old,sp_v,a_b,tp_m,tind_m)
    
    spcid = intersect(sp_v,cid);
    spcid_ind = indexin(spcid,[1;sp_v]);
    
    #Richness measures
    rich = length(cid);
    sprich = length(spcid);
    obrich = rich - sprich;
    
    #Connectance
    conn = sum(tp_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    conn_ind = sum(tind_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    
    
    #Turnover
    turnover = 1 - (length(intersect(cid,cid_old)) / length(union(cid,cid_old)));
    
    #Resource overlap
    res_overlap = roverlap(cid,sp_v,a_b);
    #mean resource overlap
    mres_overlap = mean(res_overlap[isfinite.(res_overlap)]);
    
    return(rich,sprich,turnover,mres_overlap,res_overlap,conn,conn_ind)
    
end