# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");

reps = 1000;
S = 400;

tmax = 2000;
tswitch = 1000;

# S = 400;
probs = [
p_n=0.003,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
lambda = 0.5;

a_thresh = 0;
n_thresh = 0.2;
extmid = 0.5; #Similarity at which pr(ext) = 0.5
steep = 1.6; #higher is steeper

extinctions = [ones(Bool,tswitch);ones(Bool,tmax-tswitch)];
colonizations = [ones(Bool,tswitch);ones(Bool,tmax-tswitch)];

MaxN = convert(Int64,floor(S + S*lambda));

#Dynamic analyses
sprich = SharedArray{Int64}(reps,tmax);
lpot_col = SharedArray{Float64}(reps,tmax);
status = SharedArray{Float64}(reps,tmax);
prim_ext = SharedArray{Float64}(reps,tmax);
sec_ext = SharedArray{Float64}(reps,tmax);
cid_r = SharedArray{Bool}(reps,tmax,MaxN);

#Save a small file to record the settings of the simulation
namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/sim_settings.jld");
# namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
save(namespace,
"reps", reps,
"S", S,
"tmax", tmax,
"tswitch",tswitch,
"a_thresh", a_thresh,
"n_thresh", n_thresh,
"colonizations", colonizations,
"extinctions", extinctions,
"extmid", extmid, #Similarity at which pr(ext) = 0.5
"steep", steep,
"probs", probs,
"lambda",lambda);

@sync @parallel for r=1:reps
    
    maxsize = 0;
    
    #This will rerun the assembly process if the community does not assemble >10 species
    while maxsize < 10
        
        # int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);
        int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

        namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
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
        
        sprich[r,:],
        cid,
        lpot_col[r,:],
        status[r,:],
        prim_ext[r,:],
        sec_ext[r,:],
        cid_r[r,:,:] = assembly_trim(
            int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
            a_thresh,n_thresh,extmid,steep,colonizations,extinctions,tmax,S,MaxN);
        
        maxsize = maximum(sum(cid_r[r,:,:],2));
            
    end
        
    
    # cid_r[r,:,:] = copy(CID);
    
    #Save individually so data can be loaded in parallel
    namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    # namespace = string("$(homedir())/Dropbox/PostDoc//2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    save(namespace,
    "CID", cid_r[r,:,:]);
    
    
end
# 
# save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/sim.jld"),
# "lpot_col",lpot_col,
# "status",status,
# "prim_ext",prim_ext,
# "sec_ext",sec_ext,
# "cid_r",cid_r
# );


save(string("$(homedir())/2014_Lego/Anime/data/simbasic/sim.jld"),
"sprich",sprich,
"lpot_col",lpot_col,
"status",status,
"prim_ext",prim_ext,
"sec_ext",sec_ext,
"cid_r",cid_r
);
