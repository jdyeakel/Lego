# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");


reps = 1000;
S = 400;
tmax = 2500;
ppweight = 1/4;

a_thresh = 0.0;
n_thresh = 0.2;
extinctions = true;

probs = [
# p_n=0.04,
# p_a=0.01,
# p_m=0.04,
p_n=0.004,
p_a=0.01,
p_m=0.002,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
];


#Dynamic analyses
cid = SharedArray{}
lpot_col = SharedArray{Float64}(reps,tmax);
status = SharedArray{Float64}(reps,tmax);
prim_ext = SharedArray{Float64}(reps,tmax);
sec_ext = SharedArray{Float64}(reps,tmax);


@time @sync @parallel for r=1:reps
    
    int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);
    
    namespace = string("/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    save(namespace,
    "int_m", int_m,
    "sp_m", sp_m,
    "t_m", t_m,
    "tp_m", tp_m,
    "tind_m", tind_m,
    "mp_m", mp_m,
    "mind_m", mind_m);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    cid,
    lpot_col[r,:],
    status[r,:],
    prim_ext[r,:],
    sec_ext[r,:], = assembly_trim(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
        a_thresh,n_thresh,extinctions,tmax,S);

    spcid = intersect(sp_v,cid);
    spcid_ind = indexin(spcid,[1;sp_v]);
    
    #Trophic and degrees at tmax
    # deg,troph = structure(S,cid,sp_v,tind_m);
    
    # degrees[r,:] = vec(deg);
    # trophic[r,:] = vec(troph);
    
end
# 
# save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"),
# "rich",rich,
# "sprich",sprich,
# "turnover",turnover,
# "mres_overlap",mres_overlap,
# "conn",conn,
# "conn_ind",conn_ind,
# "prim_ext",prim_ext,
# "sec_ext",sec_ext,
# "status",status,
# "lpot_col",lpot_col,
# "avgdegree",avgdegree,
# "degrees",degrees,
# "trophic",trophic
# );


save(string("$(homedir())/2014_Lego/Anime/data/structure.jld"),
"rich",rich,
"sprich",sprich,
"turnover",turnover,
"mres_overlap",mres_overlap,
"res_overlap_dist",res_overlap_dist,
"conn",conn,
"conn_ind",conn_ind,
"prim_ext",prim_ext,
"sec_ext",sec_ext,
"status",status,
"lpot_col",lpot_col,
"avgdegree",avgdegree,
"degrees",degrees,
"trophic",trophic
);
