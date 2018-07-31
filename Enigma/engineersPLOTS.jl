@everywhere using Distributions
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD
@everywhere include("$(homedir())/2014_Lego/Enigma/src/preamble_defs.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/structure.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/dynstructure.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/trophicalc2.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/roverlap.jl")


namespace = string("$(homedir())/2014_Lego/Enigma/data/engineers/sim_settings.jld");
d1 = load(namespace);
reps = d1["reps"];
S = d1["S"];
lambdavec = d1["lambdavec"];
maxits = d1["maxits"];


llamb = length(lambdavec);
its = llamb*reps;

# seq = [collect(2:50);100;200;1000;2000;];
# tseqmax = length(seq);
# 
# rich = SharedArray{Int64}(llamb,reps,tseqmax);
# sprich = SharedArray{Int64}(llamb,reps,tseqmax);
# turnover = SharedArray{Float64}(llamb,reps,tseqmax);
# res_overlap = SharedArray{Float64}(llamb,reps,tseqmax);
# user_overlap = SharedArray{Float64}(llamb,reps,tseqmax);
# conn = SharedArray{Float64}(llamb,reps,tseqmax);
# conn_ind = SharedArray{Float64}(llamb,reps,tseqmax);
# avgdegree = SharedArray{Float64}(llamb,reps,tseqmax);
# res_overlap_dist = SharedArray{Float64}(llamb,reps,tseqmax,S);
# user_overlap_dist = SharedArray{Float64}(llamb,reps,tseqmax,S);
# degrees = SharedArray{Int64}(llamb,reps,tseqmax,S);
# trophic = SharedArray{Float64}(llamb,reps,tseqmax,S);

# 
# @sync @parallel for i = 0:(its - 1)
# 
#     #Across lambdavec
#     a = Int64(floor(i/reps)) + 1;
#     #Across reps
#     b = mod(i,reps) + 1;
# 
#     namespace_rep = string("$(homedir())/2014_Lego/Enigma/data/engineers/int_m_",a,"_",b,".jld");
# 
#     d2 = load(namespace_rep);
#     int_m = d2["int_m"];
#     tp_m = d2["tp_m"];
#     tind_m = d2["tind_m"];
#     mp_m = d2["mp_m"];
#     mind_m = d2["mind_m"];
# 
#     namespace_cid = string("$(homedir())/2014_Lego/Enigma/data/engineers/cid_",a,"_",b,".jld");
#     d3 = load(namespace_cid);
#     CID = d3["CID"];
# 
#     a_b,
#     n_b,
#     i_b,
#     m_b,
#     n_b0,
#     sp_v,
#     int_id = preamble_defs(int_m);
# 
#     #Analysis
#     for t = 1:tseqmax
#         #construct
#         tstep = seq[t];
#         cid = find(isodd,CID[:,tstep]);
#         cid_old = find(isodd,CID[:,tstep-1]); #because we have this, seq can't start at t=1;
# 
#         rich[a,b,t], sprich[a,b,t], turnover[a,b,t], res_overlap[a,b,t], user_overlap[a,b,t], res_overlap_all, user_overlap_all, conn[a,b,t], conn_ind[a,b,t] = dynstructure(cid,cid_old,sp_v,a_b,n_b0,tp_m,tind_m);     
# 
#         res_overlap_dist[a,b,t,1:length(res_overlap_all)] = res_overlap_all;
#         #Only save species user-overlap, and not object user-overlap
#         user_overlap_dist[a,b,t,1:length(user_overlap_all)] = user_overlap_all; 
# 
#         deg,troph = structure(S,cid,sp_v,tind_m);
# 
#         degrees[a,b,t,1:length(deg)] = deg;
#         trophic[a,b,t,1:length(troph)] = troph;
#         avgdegree[a,b,t] = mean(degrees[a,b,t,1:length(deg)]);
# 
# 
#     end
# 
# end

sprich = SharedArray{Int64}(llamb,reps,maxits);
rich = SharedArray{Int64}(llamb,reps,maxits);
#Extinction casades
@sync @parallel for i = 0:(its - 1)
    
    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;
    
    # namespace_rep = string("$(homedir())/2014_Lego/Enigma/data/engineers/int_m_",a,"_",b,".jld");
    
    # d2 = load(namespace_rep);
    # int_m = d2["int_m"];
    # tp_m = d2["tp_m"];
    # tind_m = d2["tind_m"];
    # mp_m = d2["mp_m"];
    # mind_m = d2["mind_m"];
    
    namespace_cid = string("$(homedir())/2014_Lego/Enigma/data/engineers/cid_",a,"_",b,".jld");
    d3 = load(namespace_cid);
    CID = d3["CID"];
    
    # a_b,
    # n_b,
    # i_b,
    # m_b,
    # n_b0,
    # sp_v,
    # int_id = preamble_defs(int_m);
    
    for t=2:maxits
        cid = CID[:,t];
        cid_old = CID[:,t-1];
        
        sprich[a,b,t] = sum(cid[1:S]);
        rich[a,b,t] = sum(cid[:]);
    end
end

msprich = Array{Float64}(maxits-1,llamb);
mrich = Array{Float64}(maxits-1,llamb);
for t=1:maxits-1
    for l=1:llamb
        msprich[t,l] = mean(sprich[l,:,t]);
        mrich[t,l] = mean(rich[l,:,t]);
    end
end


namespace = string("$(homedir())/2014_Lego/Enigma/figures/eng/sprichengineers.pdf");
timeseq = collect(1:maxits-1);
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal = colorRampPalette(rev(brewer.pal(9,"Blues")))(100)
pdf($namespace,width=8,height=7)
image.plot(y=$lambdavec,x=$timeseq,z=$(msprich),ylab='Expected num. objects/species',xlab='Time',log='x',col=pal,useRaster=TRUE)
dev.off()
"""

namespace = string("$(homedir())/2014_Lego/Enigma/figures/eng/richengineers.pdf");
timeseq = collect(1:maxits-1);
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal = colorRampPalette(rev(brewer.pal(9,"Blues")))(100)
pdf($namespace,width=8,height=7)
image.plot(y=$lambdavec,x=$timeseq,z=$(mrich),ylab='Expected num. objects+species',xlab='Time',log='x',col=pal,useRaster=TRUE)
dev.off()
"""




namespace = string("$(homedir())/2014_Lego/Enigma/data/engineers/sim_settings.jld");
d1 = load(namespace);
reps = d1["reps"];
S = d1["S"];
lambdavec = d1["lambdavec"];
maxits = d1["maxits"];


llamb = length(lambdavec);
its = llamb*reps;

lcdf = 500;
EXTCDF = SharedArray{Int64}(its,lcdf);
extratevec = SharedArray{Float64}(its,lcdf);
engineers = SharedArray{Int64}(its,maxits);
sprich = SharedArray{Int64}(its,maxits);
rich = SharedArray{Int64}(its,maxits);

@sync @parallel for i = 0:(its - 1)    
    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;
    
    ii = i + 1;
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    namespace_rep = string("$(homedir())/2014_Lego/Enigma/data/engineers/int_m_",a,"_",b,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    namespace_cid = string("$(homedir())/2014_Lego/Enigma/data/engineers/cid_",a,"_",b,".jld");
    d3 = load(namespace_cid);
    CID = d3["CID"];
    dt = d3["clock"];
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    
    
    #Calculate CDF of extinction cascade size
    
    for t=1:maxits
        sprich[ii,t] = sum(CID[1:S,t]);
        rich[ii,t] = sum(CID[:,t]);
        #how many engineers?
        spcid = find(isodd,CID[1:S,t]);
        engineers[ii,t] = sum(sum(m_b[spcid,:],2) .> 0);
    end
    
    spdiff = diff(sprich[ii,:]);
    rdiff = diff(rich[ii,:]);
    
    colpos = find(x->x>0,spdiff);
    extpos = find(x->x<0,spdiff);
    extinctions = spdiff[extpos].*-1;
    colonizations = spdiff[colpos];
    
    extrate = extinctions ./ dt[extpos];
    colrate = colonizations ./ dt[colpos];
    
    # extratevec = collect(0:0.0001:maximum(extrate));
    extratevec[ii,:] = collect(range(0,maximum(extrate)/lcdf,lcdf));
    
    extcdf = Array{Int64}(lcdf);
    for j=1:lcdf
        extcdf[j] = length(find(x->x<extratevec[ii,j],extrate))
    end
    
    EXTCDF[ii,:] = extcdf;    
    
    # if mod(r,1) == 0
    #     println("reps =",r)
    # end
end

objects = rich .- sprich;
mobj = vec(mean(objects[:,maxits-100:maxits],2));
obsort = sortperm(mobj);
meng = vec(mean(engineers[:,maxits-100:maxits],2));
engsort = sortperm(meng);

#Convert frequencies to probabilities
EXTCDFpr = Array{Float64}(reps,lcdf);
for i=1:reps
    EXTCDFpr[i,:] = Array(EXTCDF[i,:])/maximum(Array(EXTCDF[i,:]));
end


sortalg = engsort;
namespace = string("$(homedir())/2014_Lego/Enigma/figures/eng/engcdf2.pdf");
R"""
library(RColorBrewer)
pdf($namespace,width=8,height=6)
pal = colorRampPalette(rev(brewer.pal(9,"Spectral")))($reps)
plot($(extratevec[sortalg[reps],:]),$(EXTCDFpr[sortalg[reps],:]),type='l',col=paste(pal[$reps],'60',sep=''),xlim=c(0.01,1),ylim=c(0,1),log='x',xlab='Extinction rate',ylab='Cumulative Probability')
legend(x=0.6,y=0.5,legend=sapply(seq(floor(min($meng)),ceiling(max($meng)),length.out=10),floor),col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(10),cex=0.8,pch=16,bty='n',title='Num. Eng.')
"""
for i=reps-1:-1:1
    R"""
    lines($(extratevec[sortalg[i],:]),$(EXTCDFpr[sortalg[i],:]),col=paste(pal[$i],'60',sep=''),log='x')
    """
end
R"dev.off()"


sortalg = obsort;
namespace = string("$(homedir())/2014_Lego/Enigma/figures/eng/objcdf2.pdf");
R"""
library(RColorBrewer)
pdf($namespace,width=8,height=6)
pal = colorRampPalette(rev(brewer.pal(9,"Spectral")))($reps)
plot($(extratevec[sortalg[reps],:]),$(EXTCDFpr[sortalg[reps],:]),type='l',col=paste(pal[$reps],'60',sep=''),xlim=c(0.01,1),ylim=c(0,1),log='x',xlab='Extinction rate',ylab='Cumulative Probability')
legend(x=0.6,y=0.5,legend=sapply(seq(floor(min($mobj)),ceiling(max($mobj)),length.out=10),floor),col=colorRampPalette(rev(brewer.pal(9,"Spectral")))(10),cex=0.8,pch=16,bty='n',title='Num. Obj.')
"""
for i=reps-1:-1:1
    R"""
    lines($(extratevec[sortalg[i],:]),$(EXTCDFpr[sortalg[i],:]),col=paste(pal[$i],'60',sep=''),log='x')
    """
end
R"dev.off()"


