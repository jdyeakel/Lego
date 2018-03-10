loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

reps = 100;
S = 400;
tmax = 2000;
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
prim_ext = SharedArray{Float64}(reps,tmax);
sec_ext = SharedArray{Float64}(reps,tmax);
status = SharedArray{Float64}(reps,tmax);
lpot_col = SharedArray{Float64}(reps,tmax);

#Steady state analyses
# degrees = SharedArray{Array{Int64}}(reps);
# trophic = SharedArray{Float64}(reps,S);

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
    prim_ext[r,:],
    sec_ext[r,:],
    status[r,:],
    lpot_col[r,:] = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tind_m,
        a_thresh,n_thresh,extinctions,tmax);

    spcid = intersect(sp_v,cid);
    spcid_ind = indexin(spcid,[1;sp_v]);
    
    # degrees[r,:],
    # trophic[r,:] = structure(cid,sp_v,tind_m);
    
end

save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"),
"rich",rich,
"sprich",sprich,
"turnover",turnover,
"res_overlap",res_overlap,
"conn",conn,
"prim_ext",prim_ext,
"sec_ext",sec_ext,
"status",status,
"lpot_col",lpot_col,
"degrees",degrees,
"trophic",trophic
);


d = load(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"));
#This loads the dictionary
rich = d["rich"];
sprich = d["sprich"];
turnover = d["turnover"];
res_overlap = d["res_overlap"];
conn = d["conn"];
prim_ext = d["prim_ext"];
sec_ext = d["sec_ext"];
status = d["status"];
lpot_col = d["lpot_col"];
degrees = d["degrees"];