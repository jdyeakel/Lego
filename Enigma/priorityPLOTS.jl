# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncsYOG.jl");

# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
namespace = string("$(homedir())/2014_Lego/Enigma/data/priority/sim_settings.jld");
d1 = load(namespace);
reps = d1["reps"];
S = d1["S"];
maxits = d1["maxits"];
athresh = d1["athresh"];
nthresh = d1["nthresh"];
intreps = d1["intreps"];


occm = SharedArray{Float64}(reps,intreps,S);

@sync @parallel for r=1:reps
    #Read in the interaction matrix
    
    # namespace_rep = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    namespace_rep = string("$(homedir())/2014_Lego/Enigma/data/priority/int_m",r,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    for rr = 1:intreps
        # namespace_cid = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
        namespace_cid = string("$(homedir())/2014_Lego/Enigma/data/priority/cid_",r,"_",rr,".jld");
        d3 = load(namespace_cid);
        CID = d3["CID"];
    
    

        #Analysis
        occm[r,rr,:] = vec(sum(CID[1:S,:],2)/ maxits) ;
        
    
    end

end
sporder = sortperm(vec(mean(occm[1,:,:],1)),rev=true);
namespace = string("$(homedir())/2014_Lego/Enigma/figures/priority/occurance.pdf");
R"""
pdf($namespace,width=15,height=8)
boxplot($(occm[1,:,sporder]),xlab='Species ordered by persistance',ylab='Persistance')
dev.off()
"""
