# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");



#Establish community template
S = 400;
ppweight = 1/4;
# S = 400;
probs = [
# p_n=0.04,
# p_a=0.02,
# p_m=0.01,
p_n=0.1,
p_a=0.1,
# p_m=0.01,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

pevec = collect(0:0.1:5);
n = zeros(length(pevec));

nev=zeros(length(pevec))
for i=1:length(pevec)
  ne=0;
  p_engineer = pevec[i];
  int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv2(S,p_engineer,probs,ppweight);
  for j=1:S
    if length(find(x->x=='m',int_m[j,:])) > 0
      ne += 1;
    end
  end  
  n[i] = length(find(x->x=='i',diag(int_m)));
  nev[i] = ne;
end

R"plot($n)"








pm_vec = collect(0.0:0.0002:0.01);
lpm_vec = length(pm_vec);

reps = 25;
S = 400;
tmax = 2500;
ppweight = 1/4;

a_thresh = 0.0;
n_thresh = 0.2;
extinctions = true;

#Dynamic analyses
cid = SharedArray{}
lpot_col = SharedArray{Float64}(lpm_vec,reps,tmax);
status = SharedArray{Float64}(lpm_vec,reps,tmax);
prim_ext = SharedArray{Float64}(lpm_vec,reps,tmax);
sec_ext = SharedArray{Float64}(lpm_vec,reps,tmax);
cid_r = SharedArray{Bool}(lpm_vec,reps,tmax,S*2);

#Save a small file to record the settings of the simulation
namespace = string("$(homedir())/2014_Lego/Anime/data/simengineer/sim_settings.jld");
# namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
save(namespace,
"reps", reps,
"S", S,
"tmax", tmax,
"a_thresh", a_thresh,
"n_thresh", n_thresh,
"extinctions", extinctions,
"pm_vec", pm_vec);

@time @sync @parallel for r=1:reps
    
        for i=1:lpm_vec
        
        probs = [
        # p_n=0.04,
        # p_a=0.01,
        # p_m=0.04,
        p_n=0.004,
        p_a=0.01,
        p_m=pm_vec[i]#0.002,
        p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
        ];

        int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);
        
        namespace = string("$(homedir())/2014_Lego/Anime/data/simengineer/int_m_r",r,"_i",i,".jld");
        # namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
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
        lpot_col[i,r,:],
        status[i,r,:],
        prim_ext[i,r,:],
        sec_ext[i,r,:],
        CID = assembly_trim(
            int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
            a_thresh,n_thresh,extinctions,tmax,S);
        
        cid_r[i,r,:,:] = copy(CID);
        
        #Save individually so data can be loaded in parallel
        namespace = string("$(homedir())/2014_Lego/Anime/data/simengineer/cid_r",r,"_i",i,".jld");
        # namespace = string("$(homedir())/Dropbox/PostDoc//2014_Lego/Anime/data/simbasic/cid_",r,".jld");
        save(namespace,
        "CID", CID);
        
    end
    
end
# 
# save(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/sim.jld"),
# "lpot_col",lpot_col,
# "status",status,
# "prim_ext",prim_ext,
# "sec_ext",sec_ext,
# "cid_r",cid_r
# );


save(string("$(homedir())/2014_Lego/Anime/data/simengineer/sim.jld"),
"lpot_col",lpot_col,
"status",status,
"prim_ext",prim_ext,
"sec_ext",sec_ext,
"cid_r",cid_r
);


