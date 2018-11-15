if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end


reps = 5000;
S = 400;

maxits = 4000;

# S = 400;
probs = (
p_n=0.001,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
);
#expected objects per species
lambda = 0.0;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));

# cidr = SharedArray{Bool}(reps,MaxN,maxits);

#Save a small file to record the settings of the simulation
namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
# namespace = string("/$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
probs_nolabel = (probs.p_n, probs.p_a);

@save namespace reps S maxits athresh nthresh lambda probs_nolabel;


@sync @distributed for r = 1:reps
    
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
    @save namespace int_m tp_m tind_m mp_m mind_m;

    sprich,rich,clock,CID = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
        athresh,nthresh,maxits);
    
    #Save individually so data can be loaded in parallel
    namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    # namespace = string("$(homedir())//2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    @save namespace CID clock;
    
    # println("reps = ",r)
end


