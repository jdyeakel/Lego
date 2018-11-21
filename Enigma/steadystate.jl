if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end


reps = 5000;
S = 400;

maxits = 4000;

SOprobs = (
p_n=0.001,
p_a=0.003
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


#expected objects per species
lambda = 0.0;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));

# cidr = SharedArray{Bool}(reps,MaxN,maxits);

#Save a small file to record the settings of the simulation
# namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
if homedir() == "/home/z840"
    namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
else
    namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
end


# probs_nolabel = (probs.p_n, probs.p_a);

@save namespace reps S maxits athresh nthresh lambda SSprobs SOprobs OOprobs;


@sync @distributed for r = 1:reps
    
    # int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    if homedir() == "/home/z840"
        namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    else
        namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    end
    
    @save namespace int_m tp_m tind_m mp_m mind_m;

    sprich,rich,clock,CID = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
        athresh,nthresh,maxits);
    
    #Save individually so data can be loaded in parallel
    if homedir() == "/home/z840"
        namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    else
        namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    end
    @save namespace CID clock;
    
    # println("reps = ",r)
end


