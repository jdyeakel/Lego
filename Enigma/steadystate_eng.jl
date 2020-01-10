if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end


reps = 5000;


cn = pi;
ce = sqrt(2);
cp = 1.0;

#How much more does a mutualism benefit relative to the penalty of a trophic link?
cn_ce = cn/ce;
#How much more does a mutualism benefit relative to the penalty of a predator?
cn_cp = cn/cp;


S = 200;
maxits = 4000;
SOprobs = (
p_n=0.002,
p_a=0.01
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


#expected objects per species
lambda = 0;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));

# cidr = SharedArray{Bool}(reps,MaxN,maxits);

#Save a small file to record the settings of the simulation
# namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
filename = "data/steadystate_eng/sim_settings.jld";
namespace = smartpath(filename);

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
    
    filename = "data/steadystate_eng/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @save namespace int_m tp_m tind_m mp_m mind_m;

    sprich,rich,clock,CID,events = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
        athresh,nthresh,maxits,cn,ce,cp);
    
    #Save individually so data can be loaded in parallel
    filename = "data/steadystate_eng/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @save namespace CID clock events;
    
    # println("reps = ",r)
end


