@everywhere using Distributions
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD
# 
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/intmatrixv3.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/preamble_defs.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly.jl")

@everywhere include("$(homedir())/2014_Lego/Enigma/src/intmatrixv3.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/preamble_defs.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/assembly.jl")


lambdavec = collect(0:0.1:2.0)
llamb = length(lambdavec);

reps = 100;
S = 400;

maxits = 5000;

# S = 400;
probs = [
p_n=0.003,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
athresh = 0;
nthresh = 0.5;


#Save a small file to record the settings of the simulation
# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/engineers/sim_settings.jld");
namespace = string("$(homedir())/2014_Lego/Enigma/data/engineers/sim_settings.jld");
save(namespace,
"reps", reps,
"S", S,
"maxits", maxits,
"athresh", athresh,
"nthresh", nthresh,
"lambdavec",lambdavec,
"probs",probs);

its = llamb*reps;
@sync @parallel for i = 0:(its - 1)
    
    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;
    
    lambda = lambdavec[a];
        
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    # namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/engineers/int_m_",a,"_",b,".jld");
    namespace = string("$(homedir())/2014_Lego/Enigma/data/engineers/int_m_",a,"_",b,".jld");
    save(namespace,
    "int_m", int_m,
    "tp_m", tp_m,
    "tind_m", tind_m,
    "mp_m", mp_m,
    "mind_m", mind_m);

    sprich,rich,clock,CID = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
        athresh,nthresh,maxits);
    
    #Save individually so data can be loaded in parallel
    # namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/engineers/cid_",r,".jld");
    namespace = string("$(homedir())/2014_Lego/Enigma/data/engineers/cid_",a,"_",b,".jld");
    save(namespace,
    "CID", CID);
end


