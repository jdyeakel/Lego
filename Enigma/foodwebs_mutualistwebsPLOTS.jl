if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

R"library(bipartite)"
R"options(warn = -1)"
R"library(igraph)"

# seq = [collect(2:50);100;200;500;1000;2000;4000];
seq = [2000;4000];
tseqmax = length(seq);
reps = 1000;
nvec = collect(0.0:0.1:2.0);
lnvec = length(nvec);

its = lnvec*reps;

# food_nest = SharedArray{Float64}(lnvec,reps,tseqmax);
# mutual_nest = SharedArray{Float64}(lnvec,reps,tseqmax);
comb_nest = SharedArray{Float64}(lnvec,reps,tseqmax);
rich = SharedArray{Float64}(lnvec,reps,tseqmax);
sdrich = SharedArray{Float64}(lnvec,reps,tseqmax);
# food_mod = SharedArray{Float64}(lnvec,reps,tseqmax);
# mutual_mod = SharedArray{Float64}(lnvec,reps,tseqmax);
comb_mod = SharedArray{Float64}(lnvec,reps,tseqmax);

@sync @distributed for i = 0:(its - 1)

    #Across lambdavec
    v = Int64(floor(i/reps)) + 1;
    #Across reps
    r = mod(i,reps) + 1;


    filename = "data/foodwebs_mutualistwebs/sim_settings.jld";
    indices = [v];
    namespace = smartpath(filename,indices);
    @load namespace S maxits athresh nthresh lambda SSprobs SOprobs OOprobs;



    filename = "data/foodwebs_mutualistwebs/int_m.jld";
    indices = [v,r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;

    filename = "data/foodwebs_mutualistwebs/cid.jld";
    indices = [v,r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;


    #Analysis
    for t = 1:tseqmax

        #construct
        tstep = seq[t];
        cid = findall(isodd,CID[:,tstep]);
        cid_old = findall(isodd,CID[:,tstep-1]); #because we have this, seq can't start at t=1;

        # food_amatrix = tp_m[[1;cid],[1;cid]];
        # mutual_amatrix = mp_m[[1;cid],[1;cid]];
        comb_amatrix = tp_m[[1;cid],[1;cid]] .+ mp_m[[1;cid],[1;cid]];

        #CALCULATE METRICS
        # rich[v,r,t] = length(cid);
        
        #look at last 100 timesteps and take means for
        #richness and sd richness
        richs = Array{Int64}(undef,100);
        for s=1:100
            richs[s] = sum(CID[:,tstep-s]);
        end
        
        rich[v,r,t] = mean(richs);
        sdrich[v,r,t] = std(richs);

        # R"""
        # options(warn = -1)
        # library(bipartite)
        # food_nestvalue=networklevel($food_amatrix,index='NODF')
        # """
        # @rget food_nestvalue;
        # if ismissing(food_nestvalue)
        #     food_nestvalue = NaN;
        # end
        # food_nest[v,r,t] = food_nestvalue;

        R"""
        options(warn = -1)
        library(bipartite)
        comb_nestvalue=networklevel($comb_amatrix,index='NODF')
        """
        @rget comb_nestvalue;
        if ismissing(comb_nestvalue)
            comb_nestvalue = NaN;
        end
        comb_nest[v,r,t] = comb_nestvalue;

        # if sum(mutual_amatrix) > 0
        #     R"""
        #     options(warn = -1)
        #     library(bipartite)
        #     mutual_nestvalue=networklevel($mutual_amatrix,index='NODF')
        #     """
        #     @rget mutual_nestvalue;
        #     if ismissing(mutual_nestvalue)
        #         mutual_nestvalue = NaN;
        #     end
        #     mutual_nest[v,r,t] = mutual_nestvalue;
        # else
        #     mutual_nest[v,r,t] = 0;
        # end
        #
        #Modularity
        R"""
        library(igraph)
        # g_food <- graph.adjacency($food_amatrix);
        # g_mutual <- graph.adjacency($mutual_amatrix);
        g_comb <- as.undirected(graph.adjacency($comb_amatrix));

        # food_modvalue <- modularity(g_food);
        # mutual_modvalue <- modularity(g_mutual);
        wtc <- cluster_walktrap(g_comb)
        comb_modvalue <- modularity(g_comb,membership(wtc));
        """
        # @rget food_modvalue;
        # @rget mutual_modvalue;
        @rget comb_modvalue;

        # food_mod[v,r,t] = food_modvalue;
        # mutual_mod[v,r,t] = mutual_modvalue;
        comb_mod[v,r,t] = comb_modvalue;
        
        
        


    end 
    
end

filename = "data/foodwebs_mutualistwebs/nestedmod.jld";
namespace = smartpath(filename);
@save namespace seq tseqmax reps nvec lnvec rich sdrich comb_nest comb_mod;


filename = "data/foodwebs_mutualistwebs/nestedmod.jld";
namespace = smartpath(filename);
@load namespace seq tseqmax reps nvec lnvec rich sdrich comb_nest comb_mod;

teval = length(seq);
nestall = vec(comb_nest[:,:,teval]);
mnest = mean(comb_nest[:,:,teval],dims=2);
mrich = vec(rich[:,:,teval]);
mstdrich = std(rich[:,:,teval],dims=2);
nvecreshaped = repeat(nvec,outer=reps);
filename = "/figures/yog/mean_nested.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$nestall,xlab='Frequency of mutualisms',ylab='NODF',pch='.',cex=1.5,col=pal[2],ylim=c(min($nestall),3))
points($nvec,$mnest,pch=16,col='black')
lines($nvec,$mnest)
dev.off()
"""


nestreps = vec(comb_nest[:,:,teval]);
richreps = vec(rich[:,:,teval]);
filename = "/figures/mean_nestedvsize.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=6)
plot(jitter($richreps),$nestreps,xlab='Species richness',ylab='NODF',pch='.')
dev.off()
"""


nestreps1 = vec(comb_nest[1,:,teval]);
richreps1 = vec(rich[1,:,teval]);
nestreps21 = vec(comb_nest[21,:,teval]);
richreps21 = vec(rich[21,:,teval]);
filename = "/figures/mean_nestedvsize2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot(jitter($richreps1),$nestreps1,xlab='Species richness',ylab='Nestedness (NODF)',pch='.',col=pal[1],ylim=c(0,3),xlim=c(100,200),cex=2)
points(jitter($richreps21),$nestreps21,pch='.',col=pal[2],cex=2)
legend(180,3,legend=c("Low pr(n)","High pr(n)"),col=pal[1:2],cex=0.8,bty='n',pch=16)
dev.off()
"""


#Modularity

teval = length(seq);
modall = vec(comb_mod[:,:,teval]);
mmod = median(comb_mod[:,:,teval],dims=2);
mrich = vec(rich[:,:,teval]);
mstdrich = std(rich[:,:,teval],dims=2);
nvecreshaped = repeat(nvec,outer=reps)
filename = "/figures/yog/mean_modularity.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$modall,xlab='Frequency of mutualisms',ylab='modularity',pch='.',cex=1.5,col=pal[3])
points($nvec,$mmod,pch=16,col='black')
lines($nvec,$mmod)
dev.off()
"""

filename = "/figures/yog/nestmod.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=6)
plot($modall,$nestall,pch='.',cex=1.5,col='black',xlab='Modularity',ylab='NODF')
dev.off()
"""
