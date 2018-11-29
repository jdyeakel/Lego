if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

reps = 100000;
S = 200;
# S = 400;
SOprobs = (
p_n=0.002,
p_a=0.01
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);

#expected objects per species
lambda = 0.0;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));

Pconn = SharedArray{Float64}(reps);
Pconn_ind = SharedArray{Float64}(reps);
Pmutconn = SharedArray{Float64}(reps);
Pmutconn_ind = SharedArray{Float64}(reps);
Pres_overlap_dist = SharedArray{Float64}(reps,S);
Puser_overlap_dist = SharedArray{Float64}(reps,S);
Pdegrees = SharedArray{Int64}(reps,S);
Ptl = SharedArray{Float64}(reps,S);


@sync @distributed for r=1:reps
    
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

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
    
    #Connectance mutualistic
    Pmutconn[r] = sum(mp_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    Pmutconn_ind[r] = sum(mind_m[spcid_ind,spcid_ind])/(length(spcid)^2);
    
    #Resource overlap
    res_overlap,user_overlap = roverlap(cid,sp_v,a_b,n_b0);
    Pres_overlap_dist[r,1:length(res_overlap)]=res_overlap;
    Puser_overlap_dist[r,1:length(user_overlap)]=user_overlap;
    
    degrees = deleteat!(vec(sum(tind_m[[1;spcid_ind],[1;spcid_ind]],dims=2)),1);
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


filename = "data/intm_structure.jld"
namespace = smartpath(filename);

@save namespace Pconn Pconn_ind Pmutconn Pmutconn_ind Pres_overlap_dist Puser_overlap_dist Pdegrees Ptl;

# save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/intm_structure.jld"),
# "Pconn",Pconn,
# "Pconn_ind",Pconn_ind,
# "Pres_overlap_dist",Pres_overlap_dist,
# "Pdegrees",Pdegrees,
# "Ptl",Ptl
# );    

