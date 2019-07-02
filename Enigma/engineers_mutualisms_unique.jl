if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

reps = 50;
S = 200;
maxits = 4000;

nvec = collect(0.0:0.1:2.0);
lnvec = length(nvec);

lambdavec = collect(0:0.1:2.0);
llamb = length(lambdavec);

cn = pi;
ce = sqrt(2);
cp = 1.0;

lvec = copy(lnvec);
paramvec = Array{Int64}(undef,lvec*lvec*reps,3);
paramvec[:,1] = repeat(collect(1:lvec),inner=lvec*reps);
paramvec[:,2] = repeat(collect(1:lvec),inner=reps,outer=lvec);
paramvec[:,3] = repeat(collect(1:reps),outer=lvec*lvec);


@sync @distributed for ii=1:(lvec*lvec*reps)
    
    v = paramvec[ii,1];
    w = paramvec[ii,2];
    r = paramvec[ii,3];
        
    
    avalue = 0.01;
    nvalue = (avalue/10)*nvec[v];    
    
    SOprobs = (
    p_n=nvalue,
    p_a=avalue
    );
    SSmult = 1.0; OOmult = 0.0;
    SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
    OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);

    # for w=1:llamb
    #expected objects per species
    lambda = lambdavec[w];
    athresh = 0;
    nthresh = 1.0;
    MaxN = convert(Int64,floor(S + S*lambda));
    
    if r == 1
        filename = "data/engineers_mutualisms_unique/sim_settings.jld";
        indices = [v,w];
        namespace = smartpath(filename,indices);
        @save namespace reps S maxits nvec athresh nthresh lambda SSprobs SOprobs OOprobs;
    end
    
    
    # @sync @distributed for r = 1:reps
        
    # int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4_unique(S,lambda,SSprobs,SOprobs,OOprobs);
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    filename = "data/engineers_mutualisms_unique/int_m.jld";
    indices = [v,w,r];
    namespace = smartpath(filename,indices);
    @save namespace int_m tp_m tind_m mp_m mind_m;
    

    sprich,rich,clock,CID,events = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
        athresh,nthresh,maxits,cn,ce,cp);
    
    #Save individually so data can be loaded in parallel
    filename = "data/engineers_mutualisms_unique/cid.jld";
    indices = [v,w,r];
    namespace = smartpath(filename,indices);
    @save namespace CID clock events;
            
            
        # end
        
        # println(string("pn :: ",v,"/",lnvec,"; lambda :: ",w,"/",llamb));
    
    # end

end


