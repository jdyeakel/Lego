function structure()
    
    #Degree distribution
    # degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
    degrees = deleteat!(vec(sum(tind_m[[1;spcid_ind],[1;spcid_ind]],2)),1);
    #Trophic Level
    #NOTE: if you don't account for indirect-object interactions, there will be trophically disconnected species!
    # tl = trophicalc(spcid_ind,tp_m);
    tl_ind = trophicalc(spcid_ind,tind_m);
    
    conn = sum(tind_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    
    
    return(degrees,tl_ind,conn)

end