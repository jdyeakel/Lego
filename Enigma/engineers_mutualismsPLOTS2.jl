if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

filename = "data/engineers_mutualisms3/sim_settings.jld";
indices = [1,1];
namespace = smartpath(filename,indices);
@load namespace reps S maxits nvec athresh nthresh lambda SSprobs SOprobs OOprobs;

# reps = 50;
# S = 200;
# maxits = 4000;

# nvec = collect(0.0:0.1:2.0);
lnvec = length(nvec);

lambdavec = collect(0:0.05:2.0);
llamb = length(lambdavec);


numcol = SharedArray{Int64}(lnvec,llamb,reps);
numprim = SharedArray{Int64}(lnvec,llamb,reps);
numsec = SharedArray{Int64}(lnvec,llamb,reps);
numobext = SharedArray{Int64}(lnvec,llamb,reps);
clocks = SharedArray{Float64}(lnvec,llamb,reps);
mss = SharedArray{Float64}(lnvec,llamb,reps);
mprimextrate = SharedArray{Float64}(lnvec,llamb,reps);
msecextrate = SharedArray{Float64}(lnvec,llamb,reps);
mcolrate = SharedArray{Float64}(lnvec,llamb,reps);
mobextrate = SharedArray{Float64}(lnvec,llamb,reps);
mpersistance = SharedArray{Float64}(lnvec,llamb,reps);
propcol = SharedArray{Float64}(lnvec,llamb,reps);

numcol_unique = SharedArray{Int64}(lnvec,llamb,reps);
numprim_unique = SharedArray{Int64}(lnvec,llamb,reps);
numsec_unique = SharedArray{Int64}(lnvec,llamb,reps);
numobext_unique = SharedArray{Int64}(lnvec,llamb,reps);
clocks_unique = SharedArray{Float64}(lnvec,llamb,reps);
mss_unique = SharedArray{Float64}(lnvec,llamb,reps);
mprimextrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
msecextrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
mcolrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
mobextrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
mpersistance_unique = SharedArray{Float64}(lnvec,llamb,reps);
propcol_unique = SharedArray{Float64}(lnvec,llamb,reps);

lvec = copy(lnvec);
paramvec = Array{Int64}(undef,lvec*lvec*reps,3);
paramvec[:,1] = repeat(collect(1:lvec),inner=lvec*reps);
paramvec[:,2] = repeat(collect(1:lvec),inner=reps,outer=lvec);
paramvec[:,3] = repeat(collect(1:reps),outer=lvec*lvec);


@sync @distributed for ii=1:(lvec*lvec*reps)
    
    v = paramvec[ii,1];
    w = paramvec[ii,2];
    r = paramvec[ii,3];
    
    #UNIQUE
    # filename = "data/engineers_mutualisms_unique2/sim_settings.jld";
    # indices = [v,w];
    # namespace = smartpath(filename,indices);
    # @load namespace reps S maxits nvec athresh nthresh lambda SSprobs SOprobs OOprobs;
    
    # @sync @distributed for r=1:reps
        
    # filename = "data/engineers_mutualisms_unique/int_m.jld";
    # indices = [v,w,r];
    # namespace = smartpath(filename,indices);
    # @load namespace int_m tp_m tind_m mp_m mind_m;
    
    filename = "data/engineers_mutualisms_unique2/cid.jld";
    indices = [v,w,r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock events;
    CID_unique = copy(CID);
    clock_unique = copy(clock);
    events_unique = copy(events);
    
    colpos_unique = findall(x->x==1,events_unique);
    primextpos_unique = findall(x->x==1,events_unique);
    secextpos_unique = findall(x->x==2,events_unique);
    obextpos_unique = findall(x->x==3,events_unique);
    #number of primary extinctions
    numcol_unique[v,w,r] = length(colpos_unique);
    numprim_unique[v,w,r] = length(primextpos_unique);
    numsec_unique[v,w,r] = length(secextpos_unique);
    numobext_unique[v,w,r] = length(obextpos_unique);
    clocks_unique[v,w,r] = clock_unique[maxits];
    mss_unique[v,w,r] = mean(sum(CID_unique[2:S,maxits-100:maxits],dims=1));
    
    mprimextrate_unique[v,w,r] = mean(1 ./ clock_unique[primextpos_unique]);
    msecextrate_unique[v,w,r] = mean(1 ./ clock_unique[secextpos_unique]);
    mcolrate_unique[v,w,r] = mean(1 ./ clock_unique[colpos_unique]);
    mobextrate_unique[v,w,r] = mean(1 ./ clock_unique[obextpos_unique]);
    
    tstepsincom = sum(CID_unique[2:S,:],dims=2) ./ maxits;
    mpersistance_unique[v,w,r] = mean(tstepsincom[vec(findall(!iszero,tstepsincom))]);
    
    propcol_unique[v,w,r] = sum(vec(sum(CID[2:S,:],dims=2)) .> 0) / (S - 1);
    
        
    # filename = "data/engineers_mutualisms2/sim_settings.jld";
    # indices = [v,w];
    # namespace = smartpath(filename,indices);
    # @load namespace reps S maxits nvec athresh nthresh lambda SSprobs SOprobs OOprobs;
    
    # @sync @distributed for r=1:reps
        
    # filename = "data/engineers_mutualisms2/int_m.jld";
    # indices = [v,w,r];
    # namespace = smartpath(filename,indices);
    # @load namespace int_m tp_m tind_m mp_m mind_m;
    
    filename = "data/engineers_mutualisms3/cid.jld";
    indices = [v,w,r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock events;
    
    colpos = findall(x->x==1,events);
    primextpos = findall(x->x==1,events);
    secextpos = findall(x->x==2,events);
    obextpos = findall(x->x==3,events);
    
    
    numcol[v,w,r] = length(colpos);
    numprim[v,w,r] = length(primextpos);
    numsec[v,w,r] = length(secextpos);
    numobext[v,w,r] = length(obextpos);
    clocks[v,w,r] = clock[maxits];
    mss[v,w,r] = mean(sum(CID[2:S,maxits-100:maxits],dims=1));
    
    mprimextrate[v,w,r] = mean(1 ./ clock[primextpos]);
    msecextrate[v,w,r] = mean(1 ./ clock[secextpos]);
    mcolrate[v,w,r] = mean(1 ./ clock[colpos]);
    mobextrate[v,w,r] = mean(1 ./ clock[obextpos]);
    
    tstepsincom = sum(CID[2:S,:],dims=2) ./ maxits;
    mpersistance[v,w,r] = mean(tstepsincom[vec(findall(!iszero,tstepsincom))]);
    
    propcol[v,w,r] = sum(vec(sum(CID[2:S,:],dims=2)) .> 0) / (S - 1);
    
    
end


#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)
filename = "data/engineers_mutualisms3/meanrates.jld";
namespace = smartpath(filename);
@save namespace reps lambdavec llamb nvec lnvec numcol numprim numsec numobext clocks mss mprimextrate msecextrate mcolrate mobextrate mpersistance propcol numcol_unique numprim_unique numsec_unique numobext_unique clocks_unique mss_unique mpersistance_unique mprimextrate_unique msecextrate_unique mcolrate_unique mobextrate_unique propcol_unique;


#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)
filename = "data/engineers_mutualisms3/meanrates.jld";
namespace = smartpath(filename);
@load namespace reps lambdavec llamb nvec lnvec numcol numprim numsec numobext clocks mss mprimextrate msecextrate mcolrate mobextrate mpersistance propcol numcol_unique numprim_unique numsec_unique numobext_unique clocks_unique mss_unique mpersistance_unique mprimextrate_unique msecextrate_unique mcolrate_unique mobextrate_unique propcol_unique;

avalue = 0.01;
nvec_scaled = (avalue/10) .* nvec;
#mean extinction rate as a function of engineering and mutualisms
mprimextrate_surf = Array{Float64}(undef,lnvec,llamb);
msecextrate_surf = Array{Float64}(undef,lnvec,llamb);
mprimextrate_surf_unique = Array{Float64}(undef,lnvec,llamb);
msecextrate_surf_unique = Array{Float64}(undef,lnvec,llamb);
avgprimextrate_surf = Array{Float64}(undef,lnvec,llamb);
avgsecextrate_surf = Array{Float64}(undef,lnvec,llamb);
avgprimextrate_surf_unique = Array{Float64}(undef,lnvec,llamb);
avgsecextrate_surf_unique = Array{Float64}(undef,lnvec,llamb);
# mcolrate_surf = Array{Float64}(undef,lnvec,llamb);
mpersist_surf = Array{Float64}(undef,lnvec,llamb);
mpersist_surf_unique = Array{Float64}(undef,lnvec,llamb);

mss_surf = Array{Float64}(undef,lnvec,llamb);
mss_surf_unique = Array{Float64}(undef,lnvec,llamb);

propcol_surf = Array{Float64}(undef,lnvec,llamb);
propcol_surf_unique = Array{Float64}(undef,lnvec,llamb);
for v = 1:lnvec
    for w = 1:llamb
        mprimextrate_surf[v,w] = mean(mprimextrate[v,w,:]);
        msecextrate_surf[v,w] = mean(msecextrate[v,w,:]);
        mprimextrate_surf_unique[v,w] = mean(mprimextrate_unique[v,w,:]);
        msecextrate_surf_unique[v,w] = mean(msecextrate_unique[v,w,:]);
        
        avgprimextrate_surf[v,w] = mean(numprim[v,w,:] ./ clocks[v,w,:]);
        avgsecextrate_surf[v,w] = mean(numsec[v,w,:] ./ clocks[v,w,:]);
        avgprimextrate_surf_unique[v,w] = mean(numprim_unique[v,w,:] ./ clocks_unique[v,w,:]);
        avgsecextrate_surf_unique[v,w] = mean(numsec_unique[v,w,:] ./ clocks_unique[v,w,:]);
        
        
        # mcolrate_surf[v,w] = mean(mcolrate[v,w,:]);
        mpersist_surf[v,w] = mean(mpersistance[v,w,:]);
        mpersist_surf_unique[v,w] = mean(mpersistance_unique[v,w,:]);
        
        mss_surf[v,w] = mean(mss[v,w,:]);
        mss_surf_unique[v,w] = mean(mss_unique[v,w,:]);
        
        propcol_surf[v,w] = mean(propcol[v,w,:]);
        propcol_surf_unique[v,w] = mean(propcol_unique[v,w,:]);
    end
end

# 
# #x ~ row of z
# #y ~ columns of z
# filename = "figures/engmut/extrate_engmut_long.pdf";
# namespace = smartpath(filename);
# R"""
# library(RColorBrewer)
# library(fields)
# pal = brewer.pal(11,'Spectral')
# pdf($namespace,width=5,height=5)
# image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mextrate_surf)),col=pal,ylab='Frequency of mutualisms',xlab='Mean number of modifiers/species',nlevel=11)
# # contour(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mextrate_surf)), add = TRUE, drawlabels = TRUE)
# 
# dev.off()
# """
# 
# 
# #x ~ row of z
# #y ~ columns of z
# filename = "figures/engmut/persistence_engmut_long.pdf";
# namespace = smartpath(filename);
# R"""
# library(RColorBrewer)
# pal = brewer.pal(9,'YlGnBu')
# pdf($namespace,width=5,height=5)
# image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mpersist_surf)),col=pal,ylab='Frequency of mutualisms',xlab='Mean number of modifiers/species',nlevel=11)
# # contour(x=$nvec_scaled,y=$lambdavec,z=$mpersist_surf, add = TRUE, drawlabels = TRUE)
# dev.off()
# """


filename = "../manuscript/fig_engineers4.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pdf($namespace,width=8,height=6)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
   widths=c(1,1,1,1), heights=c(0.5,0.5,0.5,0.5))
par(oma = c(3, 4, 1, 1), mar = c(1, 1, 0, 1)) #,mai=c(0.6,0.6,0,0.1)
pal = brewer.pal(11,'YlOrRd')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(avgprimextrate_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Frequency of service interactions', line=2.5, cex.lab=2,las=1,outer=TRUE,xpd=NA)
title(xlab=expression(paste('Expected modifiers/species (',eta,')')), line=1.5, cex.lab=2,las=1,outer=TRUE,xpd=NA)
# title(xlab='Mean modifiers/species', line=1.5, cex.lab=1.)
mtext('1° ext.',side=4,line=0.1,las=1,padj=-10)
text(0.1,0.0019,'a',cex=2,col='Black',font=2)

pal = brewer.pal(11,'YlOrRd')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(avgsecextrate_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='Mean modifiers/species', line=1.5, cex.lab=1.)
mtext('2° ext.',side=4,line=0.1,las=1,padj=-10)
text(0.1,0.0019,'b',cex=2,col='White',font=2)

pal = brewer.pal(9,'Greens')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mpersist_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# title(ylab='', line=2.0, cex.lab=1.)
# title(xlab='', line=1.5, cex.lab=1.)
mtext('Pers.',side=4,line=0.1,las=1,padj=-10)
text(0.1,0.0019,'c',cex=2,col='Black',font=2)

pal = brewer.pal(9,'Blues')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mss_surf_unique ./ mss_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='', line=1.5, cex.lab=1.)
mtext(expression(paste({'S'['u']}^'*','/S'^'*')),side=4,line=0.1,las=1,padj=-6)
text(0.1,0.0019,'d',cex=2,col='White',font=2)

dev.off()
"""




filename = "../manuscript/fig_engineers4_unique.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pdf($namespace,width=8,height=6)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
   widths=c(1,1,1,1), heights=c(0.5,0.5,0.5,0.5))
par(oma = c(3, 4, 1, 1), mar = c(1, 1, 0, 1)) #,mai=c(0.6,0.6,0,0.1)
pal = brewer.pal(11,'YlOrRd')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(avgprimextrate_surf_unique)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Frequency of service interactions', line=2.5, cex.lab=1.2,las=1,outer=TRUE,xpd=NA)
title(xlab=expression(paste('Expected modifiers/species (',eta,')')), line=1.0, cex.lab=1.2,las=1,outer=TRUE,xpd=NA)
# title(xlab='Mean modifiers/species', line=1.5, cex.lab=1.)
mtext('1° ext.',side=4,line=0.1,las=1,padj=-10)
text(0.1,0.0019,'a',cex=2,col='Black',font=2)

pal = brewer.pal(11,'YlOrRd')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(avgsecextrate_surf_unique)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='Mean modifiers/species', line=1.5, cex.lab=1.)
mtext('2° ext.',side=4,line=0.1,las=1,padj=-10)
text(0.1,0.0019,'b',cex=2,col='White',font=2)

pal = brewer.pal(9,'Greens')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mpersist_surf_unique)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# title(ylab='', line=2.0, cex.lab=1.)
# title(xlab='', line=1.5, cex.lab=1.)
mtext('Pers.',side=4,line=0.1,las=1,padj=-10)
text(0.1,0.0019,'c',cex=2,col='Black',font=2)

pal = brewer.pal(9,'Blues')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mss_surf_unique ./ mss_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='', line=1.5, cex.lab=1.)
mtext(expression(paste({'S'['u']}^'*','/S'^'*')),side=4,line=0.1,las=1,padj=-6)
text(0.1,0.0019,'d',cex=2,col='White',font=2)

dev.off()
"""


filename = "../manuscript/fig_steadystates2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pdf($namespace,width=8,height=7)
layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), 
   widths=c(1,1,1,1), heights=c(0.5,0.5,0.5,0.5))
par(oma = c(4, 5, 1, 1), mar = c(1, 1, 0, 3)) #,mai=c(0.6,0.6,0,0.1)
pal = brewer.pal(9,'Blues')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose( mss_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE,zlim=c(0,200))
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Frequency of service interactions', line=2.5, cex.lab=1.2,las=1,outer=TRUE,xpd=NA)
title(xlab=expression(paste('Expected modifiers/species (',eta,')')), line=1.0, cex.lab=1.2,las=1,outer=TRUE,xpd=NA)
# title(xlab='Mean modifiers/species', line=1.5, cex.lab=1.)
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='', line=1.5, cex.lab=1.)
mtext(expression(paste('S'^'*')),side=4,line=0.3,las=1,padj=-8.5)
text(0.1,0.0019,'a',cex=2,col='White',font=2)


pal = brewer.pal(9,'Blues')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mss_surf_unique)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE,zlim=c(0,200))
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='', line=1.5, cex.lab=1.)
mtext(expression(paste({'S'['u']}^'*')),side=4,line=0.3,las=1,padj=-7)
text(0.1,0.0019,'b',cex=2,col='White',font=2)


pal = brewer.pal(9,'BuPu')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose( propcol_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE,zlim=c(0,1))
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Frequency of service interactions', line=2.5, cex.lab=1.2,las=1,outer=TRUE,xpd=NA)
title(xlab=expression(paste('Expected modifiers/species (',eta,')')), line=1.0, cex.lab=1.2,las=1,outer=TRUE,xpd=NA)
# title(xlab='Mean modifiers/species', line=1.5, cex.lab=1.)
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='', line=1.5, cex.lab=1.)
mtext('Prop. col.',side=4,line=0.3,las=1,padj=-12)
text(0.1,0.0019,'c',cex=2,col='White',font=2)


pal = brewer.pal(9,'BuPu')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(propcol_surf_unique)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE,zlim=c(0,1))
# axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
# title(ylab='Freq. service interactions', line=2.0, cex.lab=1.,las=1)
# title(xlab='', line=1.5, cex.lab=1.)
mtext('Prop. col.',side=4,line=0.3,las=1,padj=-12)
text(0.1,0.0019,'d',cex=2,col='White',font=2)


dev.off()
"""


# 
# dev.off()
# """
# 
# 
# filename = "../manuscript/fig_propcol.pdf";
# namespace = smartpath(filename);
# R"""
# library(RColorBrewer)
# library(fields)
# pdf($namespace,width=8,height=4)
# layout(matrix(c(1,2), 1, 2, byrow = TRUE), 
#    widths=c(1,1), heights=c(0.5,0.5))
# par(oma = c(3, 4, 1, 1), mar = c(1, 1, 0, 1)) #,mai=c(0.6,0.6,0,0.1)