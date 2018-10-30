function structure(S,com_id,sp_v,tind_m)
    
    spcid = intersect(sp_v,com_id);
    spcid_ind = indexin(spcid,[1;sp_v]);
    
    #Degree distribution
    # degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
    degrees = deleteat!(vec(sum(tind_m[[1;spcid_ind],[1;spcid_ind]],dims=2)),1);
    #Trophic Level
    #NOTE: if you don't account for indirect-object interactions, there will be trophically disconnected species!
    # tl = trophicalc(spcid_ind,tp_m);
    tl = trophicalc(spcid_ind,tind_m);
    
    # conn = sum(tind_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    
    #Size difference of community
    remainder = zeros(S-length(tl));
    
    degrees_out = [degrees;convert(Array{Int64},remainder)];
    tl_out = [tl;convert(Array{Float64},remainder)];
    
    return(degrees_out,tl_out)

end