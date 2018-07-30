@everywhere using Distributions
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD
# 
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/intmatrixv3.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/preamble_defs.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/structure.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/trophicalc2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/roverlap.jl")

reps = 50000;
S = 400;
# S = 400;
probs = [
p_n=0.003,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
lambda = 0.5;

Pconn = SharedArray{Float64}(reps);
Pconn_ind = SharedArray{Float64}(reps);
Pres_overlap_dist = SharedArray{Float64}(reps,S);
Pdegrees = SharedArray{Int64}(reps,S);
Ptl = SharedArray{Float64}(reps,S);


@time @sync @parallel for r=1:reps
    
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    #get characteristics of community Pool
    cid = collect(2:size(int_m)[1]);
    spcid = intersect(sp_v,cid);
    spcid_ind = indexin(spcid,[1;sp_v]);
    
    #Connectance
    Pconn[r] = sum(tp_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    Pconn_ind[r] = sum(tind_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    
    #Resource overlap
    res_overlap,user_overlap = roverlap(cid,sp_v,a_b,n_b0);
    Pres_overlap_dist[r,1:length(res_overlap)]=res_overlap;
    
    degrees = deleteat!(vec(sum(tind_m[[1;spcid_ind],[1;spcid_ind]],2)),1);
    #Trophic Level
    #NOTE: if you don't account for indirect-object interactions, there will be trophically disconnected species!
    # tl = trophicalc(spcid_ind,tp_m);
    tl = trophicalc(spcid_ind,tind_m);
    
    # conn = sum(tind_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    
    #Size difference of community
    remainder = zeros(S-length(tl));
    
    Pdegrees[r,:] = [degrees;convert(Array{Int64},remainder)];
    Ptl[r,:] = [tl;convert(Array{Float64},remainder)];
    
end

save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/intm_structure.jld"),
"Pconn",Pconn,
"Pconn_ind",Pconn_ind,
"Pres_overlap_dist",Pres_overlap_dist,
"Pdegrees",Pdegrees,
"Ptl",Ptl
);    
