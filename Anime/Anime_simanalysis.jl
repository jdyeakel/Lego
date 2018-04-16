loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");
# loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");




namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/sim.jld");
# namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/sim.jld");
d1 = load(namespace);
lpot_col = d1["lpot_col"];
status = d1["status"];
prim_ext = d1["prim_ext"];
sec_ext = d1["sec_ext"];
cid_r = d1["cid_r"];

reps = size(cid_r)[1];
tmax = size(cid_r)[2];
S = convert(Int64,size(cid_r)[3]/2);

for r=1:reps
    #Read in the interaction matrix
    namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    # namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    d2 = load(namespace);
    int_m = d2["int_m"];
    sp_m = d2["sp_m"];
    t_m = d2["t_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    CID = cid_r[r,:,:];
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
        
    #Analysis
    for t = 2:tmax
        
        #construct
        cid = find(isodd,CID[t,:]);
        cid_old = find(isodd,CID[t-1,:]);
        
        
        rich[t], sprich[t], turnover[t], res_overlap[t], res_overlap_all, conn[t], conn_ind[t] = dynstructure(cid,cid_old,sp_v,a_b,tp_m,tind_m);      
        
        res_overlap_dist[t,1:length(res_overlap_all)] = res_overlap_all;
        #Trophic and degrees at tmax
        deg,troph = structure(S,cid,sp_v,tind_m);
        degrees[t,1:length(deg)] = deg;
        trophic[t,1:length(troph)] = troph;
        avgdegree[t] = mean(degrees[t,1:length(deg)]);