if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end
#git

lambdavec = collect(0:0.1:5.0);
llamb = length(lambdavec);

reps = 1000;
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
# lambda = 0.0;
athresh = 0;
nthresh = 1.0;
# MaxN = convert(Int64,floor(S + S*lambda));


#Save a small file to record the settings of the simulation
# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/engineers/sim_settings.jld");
filename = "data/engineers/sim_settings.jld";
namespace = smartpath(filename);
@save namespace reps S maxits athresh nthresh lambdavec SSprobs SOprobs OOprobs;


its = llamb*reps;
@sync @distributed for i = 0:(its - 1)
    
    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;
    
    lambda = lambdavec[a];
        
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    filename = "data/engineers/int_m.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices); 
    @save namespace int_m tp_m tind_m mp_m mind_m;
    
    sprich,rich,clock,CID = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
        athresh,nthresh,maxits);

    #Save individually so data can be loaded in parallel
    filename = "data/engineers/cid.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices);
    @save namespace CID clock;
    
end
