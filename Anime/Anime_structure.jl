loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

reps = 1000;
S = 400;
tmax = 5000;
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

#Steady state analyses
degrees = SharedArray{Int64}(reps,S);
trophic = SharedArray{Float64}(reps,S);

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
    conn[r,:],
    conn_ind[r,:],
    prim_ext[r,:],
    sec_ext[r,:],
    status[r,:],
    lpot_col[r,:],
    avgdegree[r,:] = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
        a_thresh,n_thresh,extinctions,tmax);

    spcid = intersect(sp_v,cid);
    spcid_ind = indexin(spcid,[1;sp_v]);
    
    #Trophic and degrees at tmax
    deg,troph = structure(S,cid,sp_v,tind_m);
    
    degrees[r,:] = vec(deg);
    trophic[r,:] = vec(troph);
    
end

save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"),
"rich",rich,
"sprich",sprich,
"turnover",turnover,
"mres_overlap",mres_overlap,
"conn",conn,
"conn_ind",conn_ind,
"prim_ext",prim_ext,
"sec_ext",sec_ext,
"status",status,
"lpot_col",lpot_col,
"degrees",degrees,
"trophic",trophic,
"avgdegree",avgdegree
);


d = load(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"));
#This loads the dictionary
rich = d["rich"];
sprich = d["sprich"];
turnover = d["turnover"];
mres_overlap = d["mres_overlap"];
conn = d["conn"];
prim_ext = d["prim_ext"];
sec_ext = d["sec_ext"];
status = d["status"];
lpot_col = d["lpot_col"];
# degrees = d["degrees"];
# trophic = d["trophic"];