# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");


# reps = 1000;

# S = 400;
# tmax = 2500;
# ppweight = 1/4;
# 
# a_thresh = 0.0;
# n_thresh = 0.2;
# extinctions = true;
# 
# probs = [
# # p_n=0.04,
# # p_a=0.01,
# # p_m=0.04,
# p_n=0.004,
# p_a=0.01,
# p_m=0.002,
# p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
# ];

reps = 5000;
S = 400;
ppweight = 1/4;
# S = 400;
probs = [
p_n=0.01,
p_a=0.01
]
#expected objects per species
lambda = 1/2;

a_thresh = 0.0;
n_thresh = 0.1;
tmax = 3000;
tswitch = 1500;
extinctions = [ones(Bool,tswitch);ones(Bool,tswitch)];
colonizations = [ones(Bool,tswitch);zeros(Bool,tswitch)];



#Dynamic analyses
lpot_col = SharedArray{Float64}(reps,tmax);
status = SharedArray{Float64}(reps,tmax);
prim_ext = SharedArray{Float64}(reps,tmax);
sec_ext = SharedArray{Float64}(reps,tmax);
cid_r = SharedArray{Bool}(reps,tmax,S*2);

#Save a small file to record the settings of the simulation
namespace = string("$(homedir())/2014_Lego/Anime/data/colcutoff/sim_settings.jld");
# namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
save(namespace,
"reps", reps,
"S", S,
"tmax", tmax,
"a_thresh", a_thresh,
"n_thresh", n_thresh,
"colonizations", colonizations,
"extinctions", extinctions,
"ppweight", ppweight,
"probs", probs,
"lambda",lambda,
"tswitch",tswitch);

@sync @parallel for r=1:reps
    
    # int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv2(S,lambda,probs,ppweight);
    
    namespace = string("$(homedir())/2014_Lego/Anime/data/colcutoff/int_m",r,".jld");
    # namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    save(namespace,
    "int_m", int_m,
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
    sec_ext[r,:],
    CID = assembly_trim(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
        a_thresh,n_thresh,colonizations,extinctions,tmax,S);
    
    cid_r[r,:,:] = copy(CID);
    
    #Save individually so data can be loaded in parallel
    namespace = string("$(homedir())/2014_Lego/Anime/data/colcutoff/cid_",r,".jld");
    # namespace = string("$(homedir())/Dropbox/PostDoc//2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    save(namespace,
    "CID", CID);
    
    
end
# 
# save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/sim.jld"),
# "lpot_col",lpot_col,
# "status",status,
# "prim_ext",prim_ext,
# "sec_ext",sec_ext,
# "cid_r",cid_r
# );


save(string("$(homedir())/2014_Lego/Anime/data/colcutoff/sim.jld"),
"lpot_col",lpot_col,
"status",status,
"prim_ext",prim_ext,
"sec_ext",sec_ext,
"cid_r",cid_r
);
