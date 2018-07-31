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


namespace = string("$(homedir())/2014_Lego/Enigma/figures/sprichengineers.pdf");
timeseq = collect(1:maxits-1);
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal = colorRampPalette(rev(brewer.pal(9,"Blues")))(100)
pdf($namespace,width=8,height=7)
image.plot(y=$lambdavec,x=$timeseq,z=$(msprich),ylab='Expected num. objects/species',xlab='Time',log='x',col=pal)
dev.off()
"""

namespace = string("$(homedir())/2014_Lego/Enigma/figures/richengineers.pdf");
timeseq = collect(1:maxits-1);
R"""
library(RColorBrewer)
library(fields)
pal = rev(brewer.pal(9,"Blues"))
pal = colorRampPalette(rev(brewer.pal(9,"Blues")))(100)
pdf($namespace,width=8,height=7)
image.plot(y=$lambdavec,x=$timeseq,z=$(mrich),ylab='Expected num. objects+species',xlab='Time',log='x',col=pal)
dev.off()
"""
