function trophicalc(spcid_ind,M)
    L = M[[1;spcid_ind],[1;spcid_ind]];
    R"""
    library(MASS)  
    library(NetIndices)
    rtl<-TrophInd(t($L))
    """
    @rget rtl;
    #Delete the primary resource and set the lowest trophic level to '1'
    tl = deleteat!(Array(rtl[1]),1) .- 1;
    return(tl);
end
