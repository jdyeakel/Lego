# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncs.jl");


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
rich = SharedArray{Int64}(reps,tmax);
sprich = SharedArray{Int64}(reps,tmax);
turnover = SharedArray{Float64}(reps,tmax);
mres_overlap = SharedArray{Float64}(reps,tmax);
conn = SharedArray{Float64}(reps,tmax);
conn_ind = SharedArray{Float64}(reps,tmax);
prim_ext = SharedArray{Float64}(reps,tmax);
sec_ext = SharedArray{Float64}(reps,tmax);
status = SharedArray{Float64}(reps,tmax);
lpot_col = SharedArray{Float64}(reps,tmax);
avgdegree = SharedArray{Float64}(reps,tmax);

res_overlap_dist = SharedArray{Float64}(reps,tmax,S);

#Steady state analyses
degrees = SharedArray{Int64}(reps,tmax,S);
trophic = SharedArray{Float64}(reps,tmax,S);

@time @sync @parallel for r=1:reps
    
    int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    cid,
    rich[r,:],
    sprich[r,:],
    turnover[r,:],
    mres_overlap[r,:],
    res_overlap_dist[r,:],
    conn[r,:],
    conn_ind[r,:],
    prim_ext[r,:],
    sec_ext[r,:],
    status[r,:],
    lpot_col[r,:],
    avgdegree[r,:],
    degrees[r,:,:],
    trophic[r,:,:] = assembly(
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
