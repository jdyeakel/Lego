# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");
loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");

# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/sim.jld");
namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/sim.jld");
d1 = load(namespace);
lpot_col = d1["lpot_col"];
status = d1["status"];
prim_ext = d1["prim_ext"];
sec_ext = d1["sec_ext"];
cid_r = d1["cid_r"];

reps = size(cid_r)[1];
tmax = size(cid_r)[2];
S = convert(Int64,size(cid_r)[3]/2);

#We need to 'grab' individual cid_r's within the parallel loop
#Export them here and then import them separately under parallel rep loop
for r=1:reps
    namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    # namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    CID = cid_r[r,:,:];
    save(namespace,
    "CID", CID);
end

#define cid_r everywhere


# NOTE: just analyze the timesteps necessary for the sequence

seq = [10;25;50;100;200;1000;2000;];
tseqmax = length(seq);

rich = SharedArray{Int64}(reps,tseqmax);
sprich = SharedArray{Int64}(reps,tseqmax);
turnover = SharedArray{Float64}(reps,tseqmax);
res_overlap = SharedArray{Float64}(reps,tseqmax);
conn = SharedArray{Float64}(reps,tseqmax);
conn_ind = SharedArray{Float64}(reps,tseqmax);
avgdegree = SharedArray{Float64}(reps,tseqmax);
res_overlap_dist = SharedArray{Float64}(reps,tseqmax,S);
degrees = SharedArray{Int64}(reps,tseqmax,S);
trophic = SharedArray{Float64}(reps,tseqmax,S);

@time @sync @parallel for r=1:reps
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    namespace_rep = string("$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    sp_m = d2["sp_m"];
    t_m = d2["t_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    namespace_cid = string("$(homedir())/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    d3 = load(namespace_cid);
    CID = d3["CID"];
    
    CID = copy(cid_r[r,:,:]);
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    richt = Array{Int64}(tseqmax);
    spricht = Array{Int64}(tseqmax);
    turnovert = Array{Float64}(tseqmax);
    res_overlapt = Array{Float64}(tseqmax);
    connt = Array{Float64}(tseqmax);
    conn_indt = Array{Float64}(tseqmax);
    avgdegreet = Array{Float64}(tseqmax);
    degreest = Array{Int64}(tseqmax,S);
    trophic = Array{Float64}(tseqmax,S);
    
    res_overlap_distt = Array{Float64}(tseqmax,S);

    #Analysis
    for t = 1:tseqmax
        
        #construct
        tstep = seq[t];
        cid = find(isodd,CID[tstep,:]);
        cid_old = find(isodd,CID[tstep-1,:]);
        
        richt[t], spricht[t], turnovert[t], res_overlapt[t], res_overlap_all, connt[t], conn_indt[t] = dynstructure(cid,cid_old,sp_v,a_b,tp_m,tind_m);     
        
        res_overlap_distt[t,1:length(res_overlap_all)] = res_overlap_all; 
        
        deg,troph = structure(S,cid,sp_v,tind_m);
        
        degreest[t,1:length(deg)] = deg;
        trophict[t,1:length(troph)] = troph;
        avgdegreet[t] = mean(degrees[r,t,1:length(deg)]);
        
    end
    
    rich[r,:] = richt;
    sprich[r,:] = spricht;
    turnover[r,:] = turnovert;
    res_overlap[r,:] = res_overlapt;
    
    conn[r,:] = connt;
    conn_ind[r,:] = conn_indt;
    degrees[r,:,:] = degreest;
    trophic[r,:,:] = trophict;
    avgdegree[r,:,:] = avgdegreet;
    
    res_overlap_dist[r,:,:] = res_overlap_distt;
    
end












test = SharedArray{Bool}(reps);

@time @sync @parallel for r=1:reps
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/sim.jld");
    d1 = load(namespace);
    lpot_col = d1["lpot_col"];
    status = d1["status"];
    prim_ext = d1["prim_ext"];
    sec_ext = d1["sec_ext"];
    cid_r = d1["cid_r"];
    
    namespace_rep = string("$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    sp_m = d2["sp_m"];
    t_m = d2["t_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    x = cid_r[r,1,1];
    
    test[r] = x;
end




