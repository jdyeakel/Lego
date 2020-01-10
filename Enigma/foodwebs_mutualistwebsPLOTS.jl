if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

# R"library(bipartite)"
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
mut_nest = SharedArray{Float64}(lnvec,reps,tseqmax);

# comb_nest_r = SharedArray{Float64}(lnvec,reps,tseqmax);
rich = SharedArray{Float64}(lnvec,reps,tseqmax);
sdrich = SharedArray{Float64}(lnvec,reps,tseqmax);
# food_mod = SharedArray{Float64}(lnvec,reps,tseqmax);
# mutual_mod = SharedArray{Float64}(lnvec,reps,tseqmax);
# comb_mod = SharedArray{Float64}(lnvec,reps,tseqmax);

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
        
        #trim matrix to eliminate non-mutualistic interactions
        
        
        

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
        # options(warn = -1)
        # library(bipartite)
        # comb_nestvalue=networklevel($comb_amatrix,index='NODF')
        
        # if(!require(devtools)){install.packages('devtools'); library(devtools)} 
        # install_bitbucket(“maucantor/unodf”); 
        library(UNODF)
        library(igraph)
        unodfvalue = suppressWarnings(unodf($comb_amatrix,selfloop=FALSE))
        comb_nestvalue_c = unodfvalue$UNODFc
        comb_nestvalue_r = unodfvalue$UNODFr
        """
        @rget comb_nestvalue_c;
        @rget comb_nestvalue_r;
        
        mrowsums = findall(x->x>0,vec(sum(mp_m[[1;cid],[1;cid]],dims=2)));
        mcolsums = findall(x->x>0,vec(sum(mp_m[[1;cid],[1;cid]],dims=2)));
        mutualisticspecies = unique([mrowsums;mcolsums]);
        if length(mutualisticspecies) > 0
            mut_amatrix = comb_amatrix[mutualisticspecies,mutualisticspecies];
            R"""
            mut_unodfvalue = suppressWarnings(unodf($mut_amatrix,selfloop=TRUE))
            mut_nestvalue_c = mut_unodfvalue$UNODFc
            mut_nestvalue_r = mut_unodfvalue$UNODFr
            """
            @rget mut_nestvalue_c;
            @rget mut_nestvalue_r;
        else 
            mut_nestvalue_c = 0;
            mut_nestvalue_r = 0;
        end
        
        
        
        
        # if ismissing(comb_nestvalue)
        #     comb_nestvalue = NaN;
        # end
        comb_nest[v,r,t] = mean([comb_nestvalue_c,comb_nestvalue_r]);
        mut_nest[v,r,t] = mean([mut_nestvalue_c,mut_nestvalue_r]);
        # comb_nest_r[v,r,t] = comb_nestvalue_r;

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
        # R"""
        # 
        # # g_food <- graph.adjacency($food_amatrix);
        # # g_mutual <- graph.adjacency($mutual_amatrix);
        # g_comb <- as.undirected(graph.adjacency($comb_amatrix));
        # 
        # # food_modvalue <- modularity(g_food);
        # # mutual_modvalue <- modularity(g_mutual);
        # wtc <- cluster_walktrap(g_comb)
        # comb_modvalue <- modularity(g_comb,membership(wtc));
        # """
        # # @rget food_modvalue;
        # # @rget mutual_modvalue;
        # @rget comb_modvalue;
        # 
        # # food_mod[v,r,t] = food_modvalue;
        # # mutual_mod[v,r,t] = mutual_modvalue;
        # comb_mod[v,r,t] = comb_modvalue;
        # 
        
        
    end 
    
end



filename = "data/foodwebs_mutualistwebs/nestedmod2.jld";
namespace = smartpath(filename);
@save namespace seq tseqmax reps nvec lnvec rich sdrich comb_nest mut_nest;


filename = "data/foodwebs_mutualistwebs/nestedmod2.jld";
namespace = smartpath(filename);
@load namespace seq tseqmax reps nvec lnvec rich sdrich comb_nest mut_nest;

teval = length(seq);
nestall = vec(comb_nest[:,:,teval]);
mut_nestall = vec(mut_nest[:,:,teval]);
mnest = mean(comb_nest[:,:,teval],dims=2);
mut_mnest = mean(mut_nest[:,:,teval],dims=2);
mrich = vec(rich[:,:,teval]);
mstdrich = std(rich[:,:,teval],dims=2);
nvecreshaped = repeat(nvec,outer=reps);


filename = "/figures/yog/mean_nested_rev.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$nestall,xlab='Frequency of mutualisms',ylab='UNODF',pch='.',cex=1.5,col=pal[2],ylim=c(min($nestall),0.025))
points($nvec,$mnest,pch=16,col='black')
lines($nvec,$mnest)
dev.off()
"""

filename = "/figures/yog/mean_mutnested_rev.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$mut_nestall,xlab='Frequency of mutualisms',ylab='UNODF',pch='.',cex=1.5,col=pal[2],ylim=c(0,0.1))
points($nvec,$mut_mnest,pch=16,col='black')
lines($nvec,$mut_mnest)
dev.off()
"""




#ylim=c(min($nestall),max($nestall))

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
filename = "/figures/yog/mean_nestedvsize2_rev.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot(jitter($richreps1),$nestreps1,xlab='Species richness',ylab='Nestedness (NODF)',pch='.',col=pal[1],ylim=c(0,0.025),xlim=c(100,200),cex=2)
points(jitter($richreps21),$nestreps21,pch='.',col=pal[2],cex=2)
legend(180,0.025,legend=c("Low pr(n)","High pr(n)"),col=pal[1:2],cex=0.8,bty='n',pch=16)
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







#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)


avalue = 0.01;
nvec_scaled = (avalue/10) .* nvec;



teval = length(seq);
nestall = vec(comb_nest[:,:,teval]);
mnest = mean(comb_nest[:,:,teval],dims=2);
mrich = vec(rich[:,:,teval]);
mstdrich = std(rich[:,:,teval],dims=2);
nvecreshaped = repeat(nvec,outer=reps);

filename = "data/engineers_mutualisms/meanrates.jld";
namespace = smartpath(filename);
@load namespace reps lambdavec llamb sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations mpersistance stdpersistance

mpersist_noengin = transpose(mpersistance[:,1,:]);

filename = "../manuscript/fig_nested.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,height=5,width=5)
layout(matrix(c(1,2), 2, 1, byrow = TRUE), 
   widths=c(1,1), heights=c(0.4,0.5))
par(oma = c(0.5, 1, 1, 1), mar = c(3, 4, 0, 1),mai=c(0.6,0.6,0.0,0.0))
boxplot($mpersist_noengin, names = c(0.0,rep("",9),0.001,rep("",9),0.002),ylim=c(0.5,0.9),xlab='',ylab='',axes=FALSE,pch='.',boxwex=0.6,col='gray')
axis(2,at=seq(0.5,0.9,by=0.1),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
# axis(1,at=seq(1,21,by=1),labels=FALSE,tck=-0.015,mgp=c(0.5,0.5,0)) #at=seq(1,21,by=1),labels=c("0",rep("",9),"0.001",rep("",9),"0.002")
title(ylab='Persistence', line=2.0, cex.lab=1.2)
# title(xlab='Persistence', line=2.0, cex.lab=1.2)
plot(jitter($(nvecreshaped .* (avalue/10))),$nestall,xlab='',ylab='',pch='.',cex=1.5,col=pal[2],ylim=c(min($nestall),3),axes=FALSE)
axis(2,at=seq(0,3,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
axis(1,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Nestedness (NODF)', line=2.0, cex.lab=1.2)
title(xlab='Freq. mutualisms', line=2.0, cex.lab=1.2)
points($(nvec .* (avalue/10)),$mnest,pch=16,col='black')
lines($(nvec  .* (avalue/10)),$mnest)
dev.off()
"""
