function dynstructure(cid,cid_old,sp_v,a_b,n_b0,tp_m,tind_m,int_id,athresh,nthresh)
    
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
    res_overlap, user_overlap = roverlap(cid,sp_v,a_b,n_b0);
    #mean resource overlap
    mres_overlap = mean(res_overlap[isfinite.(res_overlap)]);
    
    muser_overlap = mean(user_overlap[isfinite.(user_overlap)]);
    
    #Potential colonizers
    pc = potcol(sp_v,int_id,cid,a_b,n_b0,athresh,nthresh);
    
    return(rich,sprich,turnover,mres_overlap,muser_overlap,res_overlap,user_overlap,conn,conn_ind,pc)
    
end
