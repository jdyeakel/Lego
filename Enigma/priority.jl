# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");


reps = 100;
intreps = 100;
S = 400;

maxits = 4000;

# S = 400;
probs = [
p_n=0.001,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
lambda = 0.5;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));

# cidr = SharedArray{Bool}(reps,MaxN,maxits);

#Save a small file to record the settings of the simulation
namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
# namespace = string("/$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
save(namespace,
"reps", reps,
"S", S,
"maxits", maxits,
"athresh", athresh,
"nthresh", nthresh,
"lambda",lambda,
"probs",probs);

@sync @parallel for r = 1:reps
    
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    # namespace = string("/$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    save(namespace,
    "int_m", int_m,
    "tp_m", tp_m,
    "tind_m", tind_m,
    "mp_m", mp_m,
    "mind_m", mind_m);
    
    for rr=1:intreps
        
        sprich,rich,clock,CID = assembly(
            int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
            athresh,nthresh,maxits);
        
        #Save individually so data can be loaded in parallel
        namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/cid_",r,"_",rr,".jld");
        # namespace = string("$(homedir())//2014_Lego/Anime/data/simbasic/cid_",r,".jld");
        save(namespace,
        "CID", CID,
        "clock",clock);
    end
    
    # println("reps = ",r)
end


