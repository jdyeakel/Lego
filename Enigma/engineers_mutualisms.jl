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

for v = 1:lnvec
    
    avalue = 0.01;
    nvalue = (avalue/10)*nvec[v];    
    
    SOprobs = (
    p_n=nvalue,
    p_a=avalue
    );
    SSmult = 1.0; OOmult = 0.0;
    SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
    OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);

    for w=1:llamb
        #expected objects per species
        lambda = lambdavec[w];
        athresh = 0;
        nthresh = 1.0;
        MaxN = convert(Int64,floor(S + S*lambda));
        
        
        filename = "data/engineers_mutualisms/sim_settings.jld";
        indices = [v,w];
        namespace = smartpath(filename,indices);
        @save namespace reps S maxits nvec athresh nthresh lambda SSprobs SOprobs OOprobs;
        
        
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
            
            filename = "data/engineers_mutualisms/int_m.jld";
            indices = [v,w,r];
            namespace = smartpath(filename,indices);
            @save namespace int_m tp_m tind_m mp_m mind_m;
            

            sprich,rich,clock,CID = assembly(
                int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
                athresh,nthresh,maxits,cn,ce,cp);
            
            #Save individually so data can be loaded in parallel
            filename = "data/engineers_mutualisms/cid.jld";
            indices = [v,w,r];
            namespace = smartpath(filename,indices);
            @save namespace CID clock;
            
            
        end
        
        println(string("pn :: ",v,"/",lnvec,"; lambda :: ",w,"/",llamb));
    
    end

end


