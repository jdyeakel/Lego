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

S = 200;

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

numprim = SharedArray{Float64}(lnvec,reps);
numsec = SharedArray{Float64}(lnvec,reps);
clocks = SharedArray{Float64}(lnvec,reps);
mpersistance = SharedArray{Float64}(lnvec,reps);

speciespersist = SharedArray{Float64}(lnvec,reps,S-1);
pindegree = SharedArray{Float64}(lnvec,reps,S-1);
poutdegree = SharedArray{Float64}(lnvec,reps,S-1);
tpindegree = SharedArray{Float64}(lnvec,reps,S-1);
tpoutdegree = SharedArray{Float64}(lnvec,reps,S-1);
mpindegree = SharedArray{Float64}(lnvec,reps,S-1);
mpoutdegree = SharedArray{Float64}(lnvec,reps,S-1);

routdegree = SharedArray{Float64}(lnvec,reps,S-1);
rindegree = SharedArray{Float64}(lnvec,reps,S-1);
troutdegree = SharedArray{Float64}(lnvec,reps,S-1);
trindegree = SharedArray{Float64}(lnvec,reps,S-1);
mroutdegree = SharedArray{Float64}(lnvec,reps,S-1);
mrindegree = SharedArray{Float64}(lnvec,reps,S-1);

@sync @distributed for i = 0:(its - 1)

    #Across lambdavec
    v = Int64(floor(i/reps)) + 1;
    #Across reps
    r = mod(i,reps) + 1;


    filename = "data/foodwebs_mutualistwebs2/sim_settings.jld";
    indices = [v];
    namespace = smartpath(filename,indices);
    @load namespace S maxits athresh nthresh lambda SSprobs SOprobs OOprobs;



    filename = "data/foodwebs_mutualistwebs2/int_m.jld";
    indices = [v,r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;

    filename = "data/foodwebs_mutualistwebs2/cid.jld";
    indices = [v,r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock events;


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
    
    #Primary and Secondary extinctions
    primextpos = findall(x->x==1,events);
    secextpos = findall(x->x==2,events);
    numprim[v,r] = length(primextpos);
    numsec[v,r] = length(secextpos);
    clocks[v,r] = clock[maxits];
    
    #persistence of each species
    tstepsincom = sum(CID[2:S,:],dims=2) ./ maxits;
    
    speciespersist[v,r,1:(S-1)] = tstepsincom;
    #POTENTIAL out-degree of each species
    poutdegree[v,r,1:(S-1)] = sum(tp_m[2:S,2:S] .+ mp_m[2:S,2:S],dims=1);
    
    #POTENTIAL in-degree of each species
    pindegree[v,r,1:(S-1)] = sum(tp_m[2:S,2:S] .+ mp_m[2:S,2:S],dims=2);
    
    tpoutdegree[v,r,1:(S-1)] = sum(tp_m[2:S,2:S],dims=1);
    #POTENTIAL in-degree of each species
    tpindegree[v,r,1:(S-1)] = sum(tp_m[2:S,2:S],dims=2);
    
    mpoutdegree[v,r,1:(S-1)] = sum(mp_m[2:S,2:S],dims=1);
    #POTENTIAL in-degree of each species
    mpindegree[v,r,1:(S-1)] = sum(mp_m[2:S,2:S],dims=2);
    
    
    mpersistance[v,r] = mean(tstepsincom[vec(findall(!iszero,tstepsincom))]);
    
    #temporal analysis
    rod = zeros(Float64,S,maxits);
    rid = zeros(Float64,S,maxits);
    trod = zeros(Float64,S,maxits);
    trid = zeros(Float64,S,maxits);
    mrod = zeros(Float64,S,maxits);
    mrid = zeros(Float64,S,maxits);
    for t=1:maxits
        spcid = findall(isodd,CID[:,t]);
        rod[spcid,t] = sum(tp_m[spcid,spcid] .+ mp_m[spcid,spcid],dims=1);
        rid[spcid,t] = sum(tp_m[spcid,spcid] .+ mp_m[spcid,spcid],dims=2);
        
        trod[spcid,t] = sum(tp_m[spcid,spcid],dims=1);
        trid[spcid,t] = sum(tp_m[spcid,spcid],dims=2);
        
        mrod[spcid,t] = sum(mp_m[spcid,spcid],dims=1);
        mrid[spcid,t] = sum(mp_m[spcid,spcid],dims=2);
        
        #species not in community - set to NAN
        absentid = setdiff(collect(1:S),spcid);
        rod[absentid,t] .= NaN;
        rid[absentid,t] .= NaN;
        trod[absentid,t] .= NaN;
        trid[absentid,t] .= NaN;
        mrod[absentid,t] .= NaN;
        mrid[absentid,t] .= NaN;
    end
    # rod[rod .== 0] .= NaN;
    # rid[rod .== 0] .= NaN;
    # trod[rod .== 0] .= NaN;
    # trid[rod .== 0] .= NaN;
    # mrod[rod .== 0] .= NaN;
    # mrid[rod .== 0] .= NaN;
    
    #average across non-nan entries
    #a zero will mean its present but its lost its interaction
    routdegree[v,r,:] = vec(meanfinite(rod,2))[2:S];
    rindegree[v,r,:] = vec(meanfinite(rid,2))[2:S];
    
    troutdegree[v,r,:] = vec(meanfinite(trod,2))[2:S];
    trindegree[v,r,:] = vec(meanfinite(trid,2))[2:S];
    
    mroutdegree[v,r,:] = vec(meanfinite(mrod,2))[2:S];
    mrindegree[v,r,:] = vec(meanfinite(mrid,2))[2:S];
    
end

#with potential degrees, you can have zeros
#with realized degrees, NaNs could either represent not-present or zero degree


filename = "data/foodwebs_mutualistwebs/nestedmod3.jld";
namespace = smartpath(filename);
@save namespace seq tseqmax reps nvec lnvec rich sdrich comb_nest mut_nest numprim numsec clocks mpersistance speciespersist poutdegree pindegree tpoutdegree tpindegree mpoutdegree mpindegree routdegree rindegree troutdegree trindegree mroutdegree mrindegree;


filename = "data/foodwebs_mutualistwebs/nestedmod3.jld";
namespace = smartpath(filename);
@load namespace seq tseqmax reps nvec lnvec rich sdrich comb_nest mut_nest numprim numsec clocks mpersistance speciespersist poutdegree pindegree tpoutdegree tpindegree mpoutdegree mpindegree routdegree rindegree troutdegree trindegree mroutdegree mrindegree;

teval = length(seq);
nestall = vec(comb_nest[:,:,teval]);
mut_nestall = vec(mut_nest[:,:,teval]);
mnest = mean(comb_nest[:,:,teval],dims=2);
mut_mnest = mean(mut_nest[:,:,teval],dims=2);
mrich = vec(rich[:,:,teval]);
mstdrich = std(rich[:,:,teval],dims=2);
nvecreshaped = repeat(nvec,outer=reps);

primextrate = vec(numprim ./ clocks);
mprimext = mean(numprim ./ clocks,dims=2);
secextrate = vec(numsec ./ clocks);
msecext = mean(numsec ./ clocks,dims=2);
persistence = vec(mpersistance);
mpersist = mean(mpersistance,dims=2);

# filename = "/figures/yog/mean_nested_rev.pdf";
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/mean_mutnested.pdf";
# namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$nestall,xlab='Frequency of mutualisms',ylab='UNODF',pch='.',cex=1.5,col=pal[2],ylim=c(min($nestall),0.025))
points($nvec,$mnest,pch=16,col='black')
lines($nvec,$mnest)
dev.off()
"""

#PRIMARY EXTINCTION RATE
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/mean_mutprimext.pdf";
# namespace = smartpath(filename);
scaled_primextrate = Int64.(floor.((primextrate .- minimum(primextrate)) ./ (maximum(primextrate) .- minimum(primextrate))*100));
scaled_primextrate[vec(findall(x->x==0,scaled_primextrate))] .= 1;
R"""
library(RColorBrewer)
pal=colorRampPalette(brewer.pal(9,"YlOrRd"))(100);
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$primextrate,xlab='Frequency of mutualisms',ylab='Primary extinction rate',pch='.',cex=1.5,col=pal[$scaled_primextrate])
points($nvec,$mprimext,pch=17,col='black')
lines($nvec,$mprimext)
dev.off()
"""

#Secondary EXTINCTION RATE
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/mean_mutsecext.pdf";
# namespace = smartpath(filename);
scaled_secextrate = Int64.(floor.((secextrate .- minimum(secextrate)) ./ (maximum(secextrate) .- minimum(secextrate))*100));
scaled_secextrate[vec(findall(x->x==0,scaled_secextrate))] .= 1;
R"""
library(RColorBrewer)
pal=colorRampPalette(brewer.pal(9,"YlOrRd"))(100);
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$secextrate,xlab='Frequency of mutualisms',ylab='Secondary extinction rate',pch='.',cex=1.5,col=pal[$scaled_secextrate])
points($nvec,$msecext,pch=18,col='black')
lines($nvec,$msecext)
dev.off()
"""

#Secondary EXTINCTION RATE
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/mean_mutpersistance.pdf";
# namespace = smartpath(filename);
scaled_persistence = Int64.(floor.((persistence .- minimum(persistence)) ./ (maximum(persistence) .- minimum(persistence))*100));
scaled_persistence[vec(findall(x->x==0,scaled_persistence))] .= 1;
R"""
library(RColorBrewer)
pal=colorRampPalette(brewer.pal(9,"Greens"))(100);
pdf($namespace,height=5,width=6)
plot(jitter($nvecreshaped),$persistence,xlab='Frequency of mutualisms',ylab='Secondary extinction rate',pch='.',cex=1.5,col=pal[$scaled_persistence])
points($nvec,$mpersist,pch=15,col='black')
lines($nvec,$mpersist)
dev.off()
"""


#all together
# filename = "/figures/yog/mean_nested_rev.pdf";


namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/mean_alltogether.pdf";
scaled_nestall = Int64.(floor.((nestall .- minimum(nestall)) ./ (maximum(nestall) .- minimum(nestall))*100));
scaled_nestall[vec(findall(x->x==0,scaled_nestall))] .= 1;
scaled_primextrate = Int64.(floor.((primextrate .- minimum(primextrate)) ./ (maximum(primextrate) .- minimum(primextrate))*100));
scaled_primextrate[vec(findall(x->x==0,scaled_primextrate))] .= 1;
scaled_secextrate = Int64.(floor.((secextrate .- minimum(secextrate)) ./ (maximum(secextrate) .- minimum(secextrate))*100));
scaled_secextrate[vec(findall(x->x==0,scaled_secextrate))] .= 1;
scaled_persistence = Int64.(floor.((persistence .- minimum(persistence)) ./ (maximum(persistence) .- minimum(persistence))*100));
scaled_persistence[vec(findall(x->x==0,scaled_persistence))] .= 1;
# namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(png)
library(fields)
pal=colorRampPalette(brewer.pal(9,"Blues")[4:9])(100);
pdf($namespace,height=8,width=4)
# layout(matrix(c(1,2,3,4),4,1,byrow=TRUE), widths=c(1,1,1,1), heights=c(0.5,0.5,0.5,0.5))
par(list(oma = c(3, 3, 0, 0), mar = c(2, 2, 1, 0)))

par(list(new=TRUE, plt=c(.075, 1, .75, 1)))
plot(jitter($nvecreshaped),$nestall * 100,xlab='',ylab='',pch='.',cex=1.5,axes=FALSE,xaxt='n',yaxt='n',col=pal[$scaled_nestall])
points($nvec,$mnest * 100,pch=16,col='black')
lines($nvec,$mnest * 100)
axis(side=2,at=seq(0,2.5,by=0.5),line=0,las=1)
mtext(side=2,'Nestedness (UNODF)',line=2.5,adj=-0.4)
# mtext(side=3,'a',font=2,outer=TRUE)

par(list(new=TRUE, plt=c(.075, 1, .48, .77)))
pal=colorRampPalette(brewer.pal(9,"YlOrRd"))(100);
plot(jitter($nvecreshaped),$primextrate,xlab='',ylab='',pch='.',cex=1.5,col=pal[$scaled_primextrate],axes=FALSE,xaxt='n',yaxt='n')
points($nvec,$mprimext,pch=17,col='black')
lines($nvec,$mprimext)
axis(side=2,at=seq(10,25,by=5),line=0,las=1)
mtext(side=2,'1° extinction rate',line=2.5)

par(list(new=TRUE, plt=c(.075, 1, .25, .5)))
pal=colorRampPalette(brewer.pal(9,"YlOrRd"))(100);
plot(jitter($nvecreshaped),$secextrate,xlab='',ylab='',pch='.',cex=1.5,col=pal[$scaled_secextrate],axes=FALSE,xaxt='n',yaxt='n')
points($nvec,$msecext,pch=18,col='black')
lines($nvec,$msecext)
axis(side=2,at=seq(2,12,by=2),line=0,las=1)
mtext(side=2,'2° extinction rate',line=2.5)


par(list(new=TRUE, plt=c(.075, 1,0, .25)))
pal=colorRampPalette(brewer.pal(9,"Greens"))(100);
plot(jitter($nvecreshaped),$persistence,xlab='',ylab='',pch='.',cex=1.5,col=pal[$scaled_persistence],axes=FALSE)
points($nvec,$mpersist,pch=15,col='black')
lines($nvec,$mpersist)
axis(side=2,at=seq(0.5,1,by=0.2),line=0,las=1)
mtext(side=2,'Persistence',line=2.5)

axis(side=1,at=seq(0,2,by=0.5),line=-1,las=1)
mtext(side=1,'Frequency of mutualisms',line=1.5)


dev.off()
"""


persistarray = meanfinite(speciespersist,3)[1:21,1:1000];

pindegreearray = meanfinite(pindegree,3)[1:21,1:1000];
poutdegreearray = meanfinite(poutdegree,3)[1:21,1:1000];

tpindegreearray = meanfinite(tpindegree,3)[1:21,1:1000];
tpoutdegreearray = meanfinite(tpoutdegree,3)[1:21,1:1000];

mpindegreearray = meanfinite(mpindegree,3)[1:21,1:1000];
mpoutdegreearray = meanfinite(mpoutdegree,3)[1:21,1:1000];


rindegreearray = meanfinite(rindegree,3)[1:21,1:1000];
routdegreearray = meanfinite(routdegree,3)[1:21,1:1000];

trindegreearray = meanfinite(trindegree,3)[1:21,1:1000];
troutdegreearray = meanfinite(troutdegree,3)[1:21,1:1000];

mrindegreearray = meanfinite(mrindegree,3)[1:21,1:1000];
mroutdegreearray = meanfinite(mroutdegree,3)[1:21,1:1000];



namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree2.pdf";
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=12,height=10)
layout(matrix(c(1,2,3,4,5,6,7,8,9),3,3,byrow=TRUE))
par(list(oma = c(3, 3, 1, 1), mar = c(5, 4, 1, 1)))

plot($(vec(tpindegreearray[1,:])),$(vec(persistarray[1,:])),pch=16,xlab='Mean trophic degree (potential)',ylab='Mean persistence - all species',col=paste(pal[1],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Potential',cex=1.5)
# text(par("usr")[1],par("usr")[3]*1.12,'Service interaction freq = 0.000',cex=1.5)
plot($(vec(tpindegreearray[10,:])),$(vec(persistarray[10,:])),pch=16,xlab='Mean trophic degree (potential)',ylab='Mean persistence - all species',col=paste(pal[1],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Potential',cex=1.5)
# text(par("usr")[1],par("usr")[3]*1.12,'Service interaction freq = 0.001',cex=1.5)
plot($(vec(tpindegreearray[21,:])),$(vec(persistarray[21,:])),pch=16,xlab='Mean trophic degree (potential)',ylab='Mean persistence - all species',col=paste(pal[1],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Potential',cex=1.5)
# text(par("usr")[1],par("usr")[3]*1.12,'Service interaction freq = 0.002',cex=1.5)


plot($(vec(trindegreearray[1,:])),$(vec(persistarray[1,:])),pch=16,xlab='Mean trophic degree (realized)',ylab='Mean persistence - colonizers',col=paste(pal[1],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Realized',cex=1.5)
# text(par("usr")[1]*1.14,par("usr")[3]*1.12,'Service interaction freq = 0.000',cex=1.5)
plot($(vec(trindegreearray[10,:])),$(vec(persistarray[10,:])),pch=16,xlab='Mean trophic degree (realized)',ylab='Mean persistence - colonizers',col=paste(pal[1],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Realized',cex=1.5)
# text(par("usr")[1]*1.14,par("usr")[3]*1.12,'Service interaction freq = 0.001',cex=1.5)
plot($(vec(trindegreearray[21,:])),$(vec(persistarray[21,:])),pch=16,xlab='Mean trophic degree (realized)',ylab='Mean persistence - colonizers',col=paste(pal[1],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Realized',cex=1.5)
# text(par("usr")[1]*1.14,par("usr")[3]*1.12,'Service interaction freq = 0.002',cex=1.5)

plot($(vec(mpindegreearray[1,:])),$(vec(persistarray[1,:])),pch=16,xlab='Mean mutualistic degree',ylab='Mean persistence - all species',col=paste(pal[2],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Potential/Realized',cex=1.5)
# text(par("usr")[1]*1.14,par("usr")[3]*1.12,'Service interaction freq = 0.000',cex=1.5)
plot($(vec(mpindegreearray[10,:])),$(vec(persistarray[10,:])),pch=16,xlab='Mean mutualistic degree',ylab='Mean persistence - all species',col=paste(pal[2],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Potential/Realized',cex=1.5)
# text(par("usr")[1]*1.14,par("usr")[3]*1.12,'Service interaction freq = 0.001',cex=1.5)
plot($(vec(mpindegreearray[21,:])),$(vec(persistarray[21,:])),pch=16,xlab='Mean mutualistic degree',ylab='Mean persistence - all species',col=paste(pal[2],'25',sep=''),ylim=c(0.3,1),cex.lab=1.5,cex=2)
# text(par("usr")[1]*1.05,par("usr")[4]*0.95,'Potential/Realized',cex=1.5)
# text(par("usr")[1]*1.14,par("usr")[3]*1.12,'Service interaction freq = 0.002',cex=1.5)

# plot($(vec(mrindegreearray[1,:])),$(vec(persistarray[1,:])),pch='.',xlab='Mean M deg. (time-avg realzed) - no needs',ylab='Mean persistence - colonizers',col=pal[2],ylim=c(0.3,1),cex.lab=1.5)
# plot($(vec(mrindegreearray[10,:])),$(vec(persistarray[10,:])),pch='.',xlab='Mean M deg. (time-avg realzed) - mid needs',ylab='Mean persistence - colonizers',col=pal[2],ylim=c(0.3,1),cex.lab=1.5)
# plot($(vec(mrindegreearray[21,:])),$(vec(persistarray[21,:])),pch='.',xlab='Mean M deg. (time-avg realzed) - max needs',ylab='Mean persistence - colonizers',col=pal[2],ylim=c(0.3,1),cex.lab=1.5)
dev.off()
"""


#boxplots for potential m-degree and persistence

namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/manuscript/fig_persistdegree_boxall.pdf";
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=10,height=6)
layout(matrix(c(1,3,2,4),2,2,byrow=TRUE))
par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
"""

m1pos_in = findall(x-> x==1,vec(trindegree));
m2pos_in = findall(x->x==2,vec(trindegree));
m3pos_in = findall(x->x==3,vec(trindegree));
m4pos_in = findall(x->x==4,vec(trindegree));
m5pos_in = findall(x->x==5,vec(trindegree));
m6pos_in = findall(x->x==6,vec(trindegree));
m7pos_in = findall(x->x==7,vec(trindegree));
m8pos_in = findall(x->x==8,vec(trindegree));
m1pos_out = findall(x-> x==1,vec(troutdegree));
m2pos_out = findall(x->x==2,vec(troutdegree));
m3pos_out = findall(x->x==3,vec(troutdegree));
m4pos_out = findall(x->x==4,vec(troutdegree));
m5pos_out = findall(x->x==5,vec(troutdegree));
m6pos_out = findall(x->x==6,vec(troutdegree));
m7pos_out = findall(x->x==7,vec(troutdegree));
m8pos_out = findall(x->x==8,vec(troutdegree));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[m1pos_in]),$(speciespersist[m2pos_in]),$(speciespersist[m3pos_in]),$(speciespersist[m4pos_in]),$(speciespersist[m5pos_in]))
# ,$(speciespersist[m6pos_in]),$(speciespersist[m7pos_in]),$(speciespersist[m8pos_in])
f <-c(rep(1,length($(speciespersist[m1pos_in]))),rep(2,length($(speciespersist[m2pos_in]))),rep(3,length($(speciespersist[m3pos_in]))),rep(4,length($(speciespersist[m4pos_in]))),rep(5,length($(speciespersist[m5pos_in]))))
#,rep(6,length($(speciespersist[m6pos_in]))),rep(7,length($(speciespersist[m7pos_in]))),rep(8,length($(speciespersist[m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[m1pos_out]),$(speciespersist[m2pos_out]),$(speciespersist[m3pos_out]),$(speciespersist[m4pos_out]),$(speciespersist[m5pos_out]))
# ,$(speciespersist[m6pos_out]),$(speciespersist[m7pos_out]),$(speciespersist[m8pos_out])
f <-c(rep(1,length($(speciespersist[m1pos_out]))),rep(2,length($(speciespersist[m2pos_out]))),rep(3,length($(speciespersist[m3pos_out]))),rep(4,length($(speciespersist[m4pos_out]))),rep(5,length($(speciespersist[m5pos_out]))))
# ,rep(6,length($(speciespersist[m6pos_out]))),rep(7,length($(speciespersist[m7pos_out]))),rep(8,length($(speciespersist[m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
# dev.off()
"""


m1pos_in = findall(x-> x==1,vec(mpindegree));
m2pos_in = findall(x->x==2,vec(mpindegree));
m3pos_in = findall(x->x==3,vec(mpindegree));
m4pos_in = findall(x->x==4,vec(mpindegree));
m5pos_in = findall(x->x==5,vec(mpindegree));
m1pos_out = findall(x-> x==1,vec(mpoutdegree));
m2pos_out = findall(x->x==2,vec(mpoutdegree));
m3pos_out = findall(x->x==3,vec(mpoutdegree));
m4pos_out = findall(x->x==4,vec(mpoutdegree));
m5pos_out = findall(x->x==5,vec(mpoutdegree));

R"""
data = c($(speciespersist[m1pos_in]),$(speciespersist[m2pos_in]),$(speciespersist[m3pos_in]),$(speciespersist[m4pos_in]),$(speciespersist[m5pos_in]))
f <-c(rep(1,length($(speciespersist[m1pos_in]))),rep(2,length($(speciespersist[m2pos_in]))),rep(3,length($(speciespersist[m3pos_in]))),rep(4,length($(speciespersist[m4pos_in]))),rep(5,length($(speciespersist[m5pos_in]))))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[2],boxwex=0.25,xlab='Mutualism in-degree (service receivers)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[m1pos_out]),$(speciespersist[m2pos_out]),$(speciespersist[m3pos_out]),$(speciespersist[m4pos_out]),$(speciespersist[m5pos_out]))
f <-c(rep(1,length($(speciespersist[m1pos_out]))),rep(2,length($(speciespersist[m2pos_out]))),rep(3,length($(speciespersist[m3pos_out]))),rep(4,length($(speciespersist[m4pos_out]))),rep(5,length($(speciespersist[m5pos_out]))))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[2],boxwex=0.25,xlab='Mutualism out-degree (service donors)',ylab='Persistence',ylim=c(0,1),las=1)
dev.off()
"""


#Trophic without mutualisms

namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/manuscript/fig_persistdegree_boxtrophic.pdf";
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=10,height=6)
layout(matrix(c(1,3,2,4),2,2,byrow=TRUE))
par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
"""

m1pos_in = findall(x-> x==1,vec(trindegree[1,:,:]));
m2pos_in = findall(x->x==2,vec(trindegree[1,:,:]));
m3pos_in = findall(x->x==3,vec(trindegree[1,:,:]));
m4pos_in = findall(x->x==4,vec(trindegree[1,:,:]));
m5pos_in = findall(x->x==5,vec(trindegree[1,:,:]));
m6pos_in = findall(x->x==6,vec(trindegree[1,:,:]));
m7pos_in = findall(x->x==7,vec(trindegree[1,:,:]));
m8pos_in = findall(x->x==8,vec(trindegree[1,:,:]));
m1pos_out = findall(x-> x==1,vec(troutdegree[1,:,:]));
m2pos_out = findall(x->x==2,vec(troutdegree[1,:,:]));
m3pos_out = findall(x->x==3,vec(troutdegree[1,:,:]));
m4pos_out = findall(x->x==4,vec(troutdegree[1,:,:]));
m5pos_out = findall(x->x==5,vec(troutdegree[1,:,:]));
m6pos_out = findall(x->x==6,vec(troutdegree[1,:,:]));
m7pos_out = findall(x->x==7,vec(troutdegree[1,:,:]));
m8pos_out = findall(x->x==8,vec(troutdegree[1,:,:]));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[m1pos_in]),$(speciespersist[m2pos_in]),$(speciespersist[m3pos_in]),$(speciespersist[m4pos_in]))
# ,$(speciespersist[m5pos_in]),$(speciespersist[m6pos_in]),$(speciespersist[m7pos_in]),$(speciespersist[m8pos_in])
f <-c(rep(1,length($(speciespersist[m1pos_in]))),rep(2,length($(speciespersist[m2pos_in]))),rep(3,length($(speciespersist[m3pos_in]))),rep(4,length($(speciespersist[m4pos_in]))))
# ,rep(5,length($(speciespersist[m5pos_in]))),rep(6,length($(speciespersist[m6pos_in]))),rep(7,length($(speciespersist[m7pos_in]))),rep(8,length($(speciespersist[m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[m1pos_out]),$(speciespersist[m2pos_out]),$(speciespersist[m3pos_out]),$(speciespersist[m4pos_out]))
# ,$(speciespersist[m5pos_out]),$(speciespersist[m6pos_out]),$(speciespersist[m7pos_out]),$(speciespersist[m8pos_out])
f <-c(rep(1,length($(speciespersist[m1pos_out]))),rep(2,length($(speciespersist[m2pos_out]))),rep(3,length($(speciespersist[m3pos_out]))),rep(4,length($(speciespersist[m4pos_out]))))
# ,rep(5,length($(speciespersist[m5pos_out]))),rep(6,length($(speciespersist[m6pos_out]))),rep(7,length($(speciespersist[m7pos_out]))),rep(8,length($(speciespersist[m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
# dev.off()
"""

m1pos_in = findall(x-> x==1,vec(tpindegree[1,:,:]));
m2pos_in = findall(x->x==2,vec(tpindegree[1,:,:]));
m3pos_in = findall(x->x==3,vec(tpindegree[1,:,:]));
m4pos_in = findall(x->x==4,vec(tpindegree[1,:,:]));
m5pos_in = findall(x->x==5,vec(tpindegree[1,:,:]));
m6pos_in = findall(x->x==6,vec(tpindegree[1,:,:]));
m7pos_in = findall(x->x==7,vec(tpindegree[1,:,:]));
m8pos_in = findall(x->x==8,vec(tpindegree[1,:,:]));
m9pos_in = findall(x->x==9,vec(tpindegree[1,:,:]));
m10pos_in = findall(x->x==10,vec(tpindegree[1,:,:]));
m1pos_out = findall(x-> x==1,vec(tpoutdegree[1,:,:]));
m2pos_out = findall(x->x==2,vec(tpoutdegree[1,:,:]));
m3pos_out = findall(x->x==3,vec(tpoutdegree[1,:,:]));
m4pos_out = findall(x->x==4,vec(tpoutdegree[1,:,:]));
m5pos_out = findall(x->x==5,vec(tpoutdegree[1,:,:]));
m6pos_out = findall(x->x==6,vec(tpoutdegree[1,:,:]));
m7pos_out = findall(x->x==7,vec(tpoutdegree[1,:,:]));
m8pos_out = findall(x->x==8,vec(tpoutdegree[1,:,:]));
m9pos_out = findall(x->x==9,vec(tpoutdegree[1,:,:]));
m10pos_out = findall(x->x==10,vec(tpoutdegree[1,:,:]));
# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[m1pos_in]),$(speciespersist[m2pos_in]),$(speciespersist[m3pos_in]),$(speciespersist[m4pos_in]),$(speciespersist[m5pos_in]),$(speciespersist[m6pos_in]),$(speciespersist[m7pos_in]),$(speciespersist[m8pos_in]),$(speciespersist[m9pos_in]),$(speciespersist[m10pos_in]))
# ,$(speciespersist[m5pos_in]),$(speciespersist[m6pos_in]),$(speciespersist[m7pos_in]),$(speciespersist[m8pos_in])
f <-c(rep(1,length($(speciespersist[m1pos_in]))),rep(2,length($(speciespersist[m2pos_in]))),rep(3,length($(speciespersist[m3pos_in]))),rep(4,length($(speciespersist[m4pos_in]))),rep(5,length($(speciespersist[m5pos_in]))),rep(6,length($(speciespersist[m6pos_in]))),rep(7,length($(speciespersist[m7pos_in]))),rep(8,length($(speciespersist[m8pos_in]))),rep(9,length($(speciespersist[m9pos_in]))),rep(10,length($(speciespersist[m10pos_in]))))
# ,rep(5,length($(speciespersist[m5pos_in]))),rep(6,length($(speciespersist[m6pos_in]))),rep(7,length($(speciespersist[m7pos_in]))),rep(8,length($(speciespersist[m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[m1pos_out]),$(speciespersist[m2pos_out]),$(speciespersist[m3pos_out]),$(speciespersist[m4pos_out]),$(speciespersist[m5pos_out]),$(speciespersist[m6pos_out]),$(speciespersist[m7pos_out]),$(speciespersist[m8pos_out]),$(speciespersist[m9pos_out]),$(speciespersist[m10pos_out]))
# ,$(speciespersist[m5pos_out]),$(speciespersist[m6pos_out]),$(speciespersist[m7pos_out]),$(speciespersist[m8pos_out])
f <-c(rep(1,length($(speciespersist[m1pos_out]))),rep(2,length($(speciespersist[m2pos_out]))),rep(3,length($(speciespersist[m3pos_out]))),rep(4,length($(speciespersist[m4pos_out]))),rep(5,length($(speciespersist[m5pos_out]))),rep(6,length($(speciespersist[m6pos_out]))),rep(7,length($(speciespersist[m7pos_out]))),rep(8,length($(speciespersist[m8pos_out]))),rep(9,length($(speciespersist[m9pos_out]))),rep(10,length($(speciespersist[m10pos_out]))))
# ,rep(5,length($(speciespersist[m5pos_out]))),rep(6,length($(speciespersist[m6pos_out]))),rep(7,length($(speciespersist[m7pos_out]))),rep(8,length($(speciespersist[m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
dev.off()
"""





################################






mutpos = 21;
#Service DONOR
tograb = findall(x->(isfinite(x) && x > 0.0), vec(mpoutdegree[mutpos,:,:]));
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_species.pdf";
R"""
library(RColorBrewer)
library(hexbin)
pal <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
pdf($namespace,width=10,height=5)
par(mfrow=c(1,2))
x <- $(mpoutdegree[mutpos,:,:][tograb])
y <- $(speciespersist[mutpos,:,:][tograb])
df <- data.frame(x,y)
h <- hexbin(df)
plot(h, colramp=pal)
"""
#Service RECEIVER
tograb = findall(x->(isfinite(x) && x > 0.0), vec(mpindegree[mutpos,:,:]));
R"""
x <- $(mpindegree[mutpos,:,:][tograb])
y <- $(speciespersist[mutpos,:,:][tograb])
df <- data.frame(x,y)
h <- hexbin(df)
plot(h, colramp=pal)
dev.off()
"""








# filename = "/figures/yog/mean_mutnested_rev.pdf";
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/mean_mutnested_rev2.pdf";
# namespace = smartpath(filename);
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

# filename = "data/engineers_mutualisms/meanrates.jld";
# namespace = smartpath(filename);
# @load namespace reps lambdavec llamb sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations mpersistance stdpersistance
# 
# mpersist_noengin = transpose(mpersistance[:,1,:]);

filename = "data/engineers_mutualisms3/meanrates.jld";
namespace = smartpath(filename);
@load namespace reps lambdavec llamb nvec lnvec numcol numprim numsec numobext clocks mss mprimextrate msecextrate mcolrate mobextrate mpersistance propcol numcol_unique numprim_unique numsec_unique numobext_unique clocks_unique mss_unique mpersistance_unique mprimextrate_unique msecextrate_unique mcolrate_unique mobextrate_unique propcol_unique;

mpersist_noengin = transpose(mpersistance[:,1,:]);
mprim_extinct_noengin = transpose(mprimextrate[:,1,:]);
msec_extinct_noengin = transpose(msecextrate[:,1,:]);


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
