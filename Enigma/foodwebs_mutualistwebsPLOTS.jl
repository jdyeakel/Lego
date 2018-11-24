if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

R"library(bipartite)"

seq = [collect(2:50);100;200;500;1000;2000;4000];
tseqmax = length(seq);

nest = SharedArray{Float64}()

for v=1:lnvec
    
    filename = "/data/foodwebs_mutualistwebs/sim_settings.jld";
    indices = [v];
    namespace = smartpath(filename,indices);
    @load namespace reps S maxits athresh nthresh lambda SSprobs SOprobs OOprobs;
    
    for r=1:reps
    
        filename = "/data/foodwebs_mutualistwebs/int_m.jld";
        indices = [v,r];
        namespace = smartpath(filename,indices);
        @load namespace int_m tp_m tind_m mp_m mind_m;
        
        filename = "/data/foodwebs_mutualistwebs/cid.jld";
        indices = [v,r];
        namespace = smartpath(filename,indices);
        @load namespace CID clock;
        
        
        #Analysis
        for t = 1:tseqmax
            
            #construct
            tstep = seq[t];
            cid = findall(isodd,CID[:,tstep]);
            cid_old = findall(isodd,CID[:,tstep-1]); #because we have this, seq can't start at t=1;
            
            food_amatrix = a_b[cid,cid];
            mutual_amatrix = 
        
            #CALCULATE METRICS
            R"nestvalue=networklevel($a_b,index="NODF");" @rget nestvalue;
            nest[v,r,t] = nestvalue;
            
        end
    end
end
