if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

filename = "data/steadystate/sim_settings.jld";
namespace = smartpath(filename);
@load namespace reps S maxits athresh nthresh lambda SSprobs SOprobs OOprobs;

# 
# d1 = load(namespace);
# reps = d1["reps"];
# S = d1["S"];
# maxits = d1["maxits"];
# athresh = d1["athresh"];
# nthresh = d1["nthresh"];


seq = [collect(2:50);100;200;500;1000;2000;4000];
tseqmax = length(seq);

rich = SharedArray{Int64}(reps,tseqmax);
sprich = SharedArray{Int64}(reps,tseqmax);
turnover = SharedArray{Float64}(reps,tseqmax);
res_overlap = SharedArray{Float64}(reps,tseqmax);
user_overlap = SharedArray{Float64}(reps,tseqmax);
conn = SharedArray{Float64}(reps,tseqmax);
conn_ind = SharedArray{Float64}(reps,tseqmax);
mutconn = SharedArray{Float64}(reps,tseqmax);
mutconn_ind = SharedArray{Float64}(reps,tseqmax);
avgdegree = SharedArray{Float64}(reps,tseqmax);
pc = SharedArray{Int64}(reps,tseqmax);
res_overlap_dist = SharedArray{Float64}(reps,tseqmax,S);
user_overlap_dist = SharedArray{Float64}(reps,tseqmax,S);
degrees = SharedArray{Int64}(reps,tseqmax,S);
trophic = SharedArray{Float64}(reps,tseqmax,S);

@sync @distributed for r=1:reps
    #Read in the interaction matrix
    
    # if homedir() == "/home/z840"
    #     namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    # else
    #     namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    # end
    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    
    # @load namespace_rep;
    # 
    # 
    # d2 = load(namespace_rep);
    # int_m = d2["int_m"];
    # tp_m = d2["tp_m"];
    # tind_m = d2["tind_m"];
    # mp_m = d2["mp_m"];
    # mind_m = d2["mind_m"];
    # 
    # if homedir() == "/home/z840"
    #     namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    # else
    #     namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    # end
    
    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;
    
    # d3 = load(namespace_cid);
    # CID = d3["CID"];
    # 
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    #Analysis
    for t = 1:tseqmax
        
        #construct
        tstep = seq[t];
        cid = findall(isodd,CID[:,tstep]);
        cid_old = findall(isodd,CID[:,tstep-1]); #because we have this, seq can't start at t=1;
        
        rich[r,t], sprich[r,t], turnover[r,t], res_overlap[r,t], user_overlap[r,t], res_overlap_all, user_overlap_all, conn[r,t], conn_ind[r,t], mutconn[r,t], mutconn_ind[r,t], pc_cid = dynstructure(cid,cid_old,sp_v,a_b,n_b0,tp_m,tind_m,mp_m,mind_m,int_id,athresh,nthresh);     
        
        pc[r,t] = length(pc_cid);
        
        res_overlap_dist[r,t,1:length(res_overlap_all)] = res_overlap_all;
        #Only save species user-overlap, and not object user-overlap
        user_overlap_dist[r,t,1:length(user_overlap_all)] = user_overlap_all; 
        
        deg,troph = structure(S,cid,sp_v,tind_m);
        
        degrees[r,t,1:length(deg)] = deg;
        trophic[r,t,1:length(troph)] = troph;
        avgdegree[r,t] = mean(degrees[r,t,1:length(deg)]);
        
        
    end

end


filename = "data/intm_structure.jld"
namespace = smartpath(filename);
@load namespace Pconn Pconn_ind Pmutconn Pmutconn_ind Pres_overlap_dist Puser_overlap_dist Pdegrees Ptl;

# @load string("$(homedir())/2014_Lego/Enigma/data/intm_structure.jld") 
# h = load(string("$(homedir())/2014_Lego/Enigma/data/intm_structure.jld"));
# Pconn = h["Pconn"];
# Pconn_ind = h["Pconn_ind"];
# Pres_overlap_dist = h["Pres_overlap_dist"];
# Pdegrees = h["Pdegrees"];
# Ptl = h["Ptl"];
Preps = size(Pconn)[1];



############
#Connectance
############
bins = [2;5;10;20;30;40;50;100;500;1000;2000];
conn_stitch,seq_stitch = sortassembly(conn,bins,seq);
meanconn = [mean(conn_stitch[findall(!isnan,conn_stitch[:,i]),i]) for i=1:length(bins)]

#mutualistic connectance
mutconn_stitch,seq_stitch = sortassembly(mutconn,bins,seq);
meanmutconn = [mean(mutconn_stitch[findall(!isnan,mutconn_stitch[:,i]),i]) for i=1:length(bins)]

filename = "figures/yog/conn_time.pdf"
namespace = smartpath(filename);

R"""
pdf($namespace,height=5,width=8)
boxplot($(conn_stitch),ylim=c(0.005,0.2),log='y',outline=FALSE,names=$(seq_stitch),
xlab='Time',ylab='Connectance',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(meanconn),ylim=c(0,0.03),pch=16)
lines($(meanconn),ylim=c(0,0.03),lwd=2)
lines(seq(0.001,3000),rep(mean($Pconn),3000),lty=3)
dev.off()
"""

filename = "figures/yog/mutconn_time.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=6)
boxplot($(mutconn_stitch),ylim=c(0.0001,0.005),log='y',outline=FALSE,names=$(seq_stitch),
xlab='Time',ylab='Connectance',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(meanmutconn),ylim=c(0,0.03),pch=16)
lines($(meanmutconn),ylim=c(0,0.03),lwd=2)
lines(seq(0.001,3000),rep(mean($Pmutconn),3000),lty=3)
dev.off()
"""


############
#Resource overlap
############
#Species pool
Preps = size(Pres_overlap_dist)[1];
Pmeanoverlap = Array{Float64}(undef,Preps);
for r=1:Preps
    Pmeanoverlap[r] = mean(Pres_overlap_dist[r,isnan.(Pres_overlap_dist[r,:]).==false]);
end

overlap_stitch,seq_stitch = sortassembly(res_overlap,bins,seq);
meanoverlap = [mean(overlap_stitch[isnan.(overlap_stitch[:,i]).==false,i]) for i=1:length(bins)];

filename = "figures/yog/trophicoverlap_time.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=6)
boxplot($(overlap_stitch),ylim=c(0,0.2),outline=FALSE,names=$(seq_stitch),
xlab='Time',ylab='Resource overlap',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(meanoverlap),ylim=c(0,0.1),pch=16)
lines($(meanoverlap),ylim=c(0,0.1),lwd=2)
lines(seq(0.001,3000),rep(mean($Pmeanoverlap),3000),lty=2)
dev.off()
"""


############
#User overlap
############
#Species pool
Preps = size(Pres_overlap_dist)[1];
Pmeanuseroverlap = Array{Float64}(undef,Preps);
for r=1:Preps
    Pmeanuseroverlap[r] = mean(Puser_overlap_dist[r,isnan.(Puser_overlap_dist[r,:]).==false]);
end

# lfseq = findall(x->x>1,diff(seq))[1];
# #This is the initial assembly process
# init_overlap = user_overlap[:,1:lfseq];
# init_overlap_trim = Array{Float64}(undef,reps,lfseq)*0;
# for r=1:reps
#     init_overlap_rm = init_overlap[r,findall(x->x>0,init_overlap[r,:])];
#     init_overlap_trim[r,1:length(init_overlap_rm)] = init_overlap_rm;
# end
# bins = [2;5;10;50;100;1000;2000;4000];
# initsteps = bins[bins.<lfseq]; #use these locations for init
# laststeps = bins[bins.>=lfseq]; #use these locations for the rest
# lastbins = indexin(laststeps,seq);

useroverlap_stitch,seq_stitch = sortassembly(user_overlap,bins,seq);
meanoverlap = [mean(useroverlap_stitch[isnan.(useroverlap_stitch[:,i]).==false,i]) for i=1:length(bins)];
filename = "figures/yog/useroverlap_time.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=6)
boxplot($(useroverlap_stitch),ylim=c(0,0.02),outline=FALSE,names=$(seq_stitch),
xlab='Time',ylab='User overlap',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(meanoverlap),ylim=c(0,0.1),pch=16)
lines($(meanoverlap),ylim=c(0,0.1),lwd=2)
#lines(seq(0.001,3000),rep(mean($Pmeanuseroverlap),3000),lty=2)
dev.off()
"""


#######################
# DEGREE AND TROPHIC DISTRIBUTION
#######################


# Somehow display how the *realized* degree distribution changes
bins = [5;10;25;50;100;200;1000;2000;4000;];
seq2 = indexin(bins,seq);

mdegt = zeros(Float64,length(seq2),S);
sddegt = zeros(Float64,length(seq2),S);
mtlt = zeros(Float64,length(seq2),S);
sdtlt = zeros(Float64,length(seq2),S);
t_tic = 0;
for t = seq2
    global t_tic = t_tic + 1;
    deg = degrees[:,t,:];
    tl = trophic[:,t,:];
    #Sort each row
    degsort = sort(deg,dims=2,rev=true);
    tlsort = sort(tl,dims=2,rev=true);
    #which column is the last non-zero?
    lastcol = findall(iszero,vec(sum(degsort,dims=1)))[1]-1;
    mdeg = Array{Float64}(undef,lastcol);
    sddeg = Array{Float64}(undef,lastcol);
    mtl = Array{Float64}(undef,lastcol);
    sdtl = Array{Float64}(undef,lastcol);
    #Take means but ignore zeros for each column through lascol
    #means and sds of degrees and trophic levels
    for i=1:lastcol
        mdeg[i] = mean(degsort[iszero.(degsort[:,i]).==false,i]);
        sddeg[i] = std(degsort[iszero.(degsort[:,i]).==false,i]);
        mtl[i] = mean(tlsort[iszero.(tlsort[:,i]).==false,i]);
        sdtl[i] = std(tlsort[iszero.(tlsort[:,i]).==false,i]);
    end
    #save the means and sds. There will be zeros for unfilled species
    mdegt[t_tic,1:length(mdeg)]=mdeg;
    sddegt[t_tic,1:length(mdeg)]=sddeg;
    mtlt[t_tic,1:length(mtl)]=mtl;
    sdtlt[t_tic,1:length(mtl)]=sdtl;
end
Pdegreesort = sort(Pdegrees,dims=2,rev=true);
Pmeandegree = vec(mapslices(mean,Pdegreesort,dims=1));
Psddeg = vec(mapslices(std,Pdegreesort,dims=1));

meanrich = round(Int64,mean(sprich[:,tseqmax]));
Pdegreesort = zeros(Int64,Preps,meanrich);
for i=1:Preps
    Pdegreesort[i,:] = sort(sample(Pdegrees[i,:],meanrich),rev=true);
end
Pmeandegree = vec(mapslices(mean,Pdegreesort,dims=1));
#If mean degree < 1, set equal to 1
firstone = findall(x->x==1,Pmeandegree)[1];
Pmeandegree[firstone:length(Pmeandegree)] .= 1;
#Where we set the degree = 1, set stdev = 0
Psddeg = vec(mapslices(std,Pdegreesort,dims=1));
Psddeg[firstone:length(Psddeg)] .= 0;

#DEGREE DISTRIBUTION

filename = "figures/degreedist_time2.pdf";
namespace = smartpath(filename);
# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/degreedist_time2.pdf");
i=length(seq2);
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
pal = brewer.pal($(length(seq2)),'Spectral')
numsp = length($(mdegt[i,iszero.(mdegt[i,:]).==false]))
plot($(mdegt[i,iszero.(mdegt[i,:]).==false]),xlim=c(1,100),ylim=c(1,6),col=pal[$i],type='l',lwd=2,xlab = 'Species sorted by degree', ylab='Mean degree')
sdev_pre = $(sddegt[i,findall(x->x>0,sddegt[i,:])]);
sdev = numeric(length($(mdegt[i,findall(x->x>0,mdegt[i,:])])))
sdev[1:length(sdev_pre)]=sdev_pre
polygon(x=c(seq(1,length(sdev)),seq(length(sdev),1)),
y=c($(mdegt[i,iszero.(mdegt[i,:]).==false])[1:length(sdev)]+sdev,
rev($(mdegt[i,iszero.(mdegt[i,:]).==false])[1:length(sdev)]-sdev)),col=paste(pal[$i],65,sep=''),border=NA)
lines($(mdegt[i,iszero.(mdegt[i,:]).==false]),xlim=c(1,200),ylim=c(0.01,50),col=pal[$i],lwd=2,xlab = 'Number of species', ylab='Median degree')
"""
for i=length(seq2)-1:-1:1
    R"""
    sdev_pre = $(sddegt[i,findall(x->x>0,sddegt[i,:])]);
    sdev = numeric(length($(mdegt[i,findall(x->x>0,mdegt[i,:])])))
    sdev[1:length(sdev_pre)]=sdev_pre
    polygon(x=c(seq(1,length(sdev)),seq(length(sdev),1)),
    y=c($(mdegt[i,iszero.(mdegt[i,:]).==false])[1:length(sdev)]+sdev,
    rev($(mdegt[i,iszero.(mdegt[i,:]).==false])[1:length(sdev)]-sdev)),col=paste(pal[$i],65,sep=''),border=NA)
    lines($(mdegt[i,iszero.(mdegt[i,:]).==false]),col=pal[$i],lwd=2)
    """
end
# R"dev.off()"

R"""
maxsp = $meanrich;
sdev = $(Psddeg)[1:maxsp];
#polygon(x=c(seq(1,maxsp),seq(maxsp,1)),
#y=c($(Pmeandegree)[1:maxsp]+sdev[1:maxsp],
#rev($(Pmeandegree)[1:maxsp]-sdev[1:maxsp])),col=paste('#000000',65,sep=''),border=NA)
lines($(Pmeandegree)[1:maxsp],col='black',type='l',lwd=2)

legend(x=215,y=120,legend = c('Pool',$(seq[seq2])),col=c('black',pal),lty=1,lwd=2,title='Time',cex=0.7,bty='n')
dev.off()
"""

#TROPHIC DISTRIBUTION

namespace = string("$(homedir())/2014_Lego/Enigma/figures/yog/trophiclevel_time.pdf");
i=length(seq2);
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
pal = brewer.pal($(length(seq2)),'Spectral')
numsp = length($(mtlt[i,!iszero.(mtlt[i,:])]))
plot($(mtlt[i,!iszero.(mtlt[i,:])]),xlim=c(1,125),ylim=c(1,10),log='y',col=pal[$i],type='l',lwd=2,xlab = 'Species sorted by trophic level', ylab='Trophic level')
sdev_pre = $(sdtlt[i,findall(x->x>0,sdtlt[i,:])]);
sdev = numeric(length($(mtlt[i,findall(x->x>0,mtlt[i,:])])))
sdev[1:length(sdev_pre)]=sdev_pre
# polygon(x=c(seq(1,length(sdev)),seq(length(sdev),1)),
# y=c($(mtlt[i,!iszero.(mtlt[i,:])])[1:length(sdev)]+sdev,
# rev($(mtlt[i,!iszero.(mtlt[i,:])])[1:length(sdev)]-sdev)),col=paste(pal[$i],65,sep=''),border=NA)
lines($(mtlt[i,!iszero.(mtlt[i,:])]),xlim=c(1,200),ylim=c(0.01,50),col=pal[$i],lwd=2,xlab = 'Number of species', ylab='Median degree')
"""
for i=length(seq2)-1:-1:1
    R"""
    sdev_pre = $(sdtlt[i,findall(x->x>0,sdtlt[i,:])]);
    sdev = numeric(length($(mtlt[i,findall(x->x>0,mtlt[i,:])])))
    sdev[1:length(sdev_pre)]=sdev_pre
    # polygon(x=c(seq(1,length(sdev)),seq(length(sdev),1)),
    # y=c($(mtlt[i,!iszero.(mtlt[i,:])])[1:length(sdev)]+sdev,
    # rev($(mtlt[i,!iszero.(mtlt[i,:])])[1:length(sdev)]-sdev)),col=paste(pal[$i],65,sep=''),border=NA)
    lines($(mtlt[i,!iszero.(mtlt[i,:])]),col=pal[$i],lwd=2)
    """
end
R"dev.off()"
# 
# 
# #####################
# # Trophic levels NOTE: TODO
# #####################
# bins = [10;25;50;100;200;1000;2000;4000;];
# seq2 = indexin(bins,seq);
# 
# R"M_t=list()"
# degsort = Array{Array{Int64}}(length(seq2));
# trophsort = Array{Array{Float64}}(length(seq2));
# for i=1:length(seq2)
#     t = seq2[i];
#     deg = reshape(degrees[:,t,:],reps*S);
#     troph = reshape(trophic[:,t,:],reps*S);
#     # deg_trim = deg[findall(!iszero,deg)];
#     # trophic_trim = troph[findall(!iszero,troph)];
#     x = vec(deg);
#     y = vec(troph);
#     zpre = [x y];
#     zkeep = findall(!iszero,sum(zpre,2));
#     z = zpre[zkeep,:];
#     sp = sortperm(z[:,1]);
#     zsort = [z[sp,1] z[sp,2]];
#     xsort = zsort[:,1];
#     ysort = zsort[:,2];
#     degsort[i] = xsort;
#     trophsort[i] = ysort;
#     R"""
#     y <- $ysort
#     x <- $xsort
#     cf <- c(0,0,0)
#     m <- try(nls(y ~ a + b * I(x^z), start = list(a = 1, b = -1, z = -1)),silent=TRUE)
#     M_t[[$i]] = m
#     """
# end
# 
# Pdeg = reshape(Array(Pdegrees),size(Pdegrees)[1]*size(Pdegrees)[2]);
# Ptroph = reshape(Array(Ptl),size(Ptl)[1]*size(Ptl)[2]);
# attached = findall(!iszero,Pdeg);
# Pdeg = Pdeg[attached];
# Ptroph = Ptroph[attached];
# z = [Pdeg Ptroph];
# sp = sortperm(z[:,1]);
# zsort = [z[sp,1] z[sp,2]];
# xsort = zsort[:,1];
# ysort = zsort[:,2];
# Pdegsort = xsort;
# Ptrophsort = ysort;
# R"""
# y <- $ysort
# x <- $xsort
# mpool <- try(nls(y ~ a + b * I(x^z), start = list(a = 1, b = -1, z = -1)),silent=TRUE)
# """
# 
# 
# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/trophic_degrees_time.pdf");
# R"""
# library(RColorBrewer)
# pdf($namespace,height=5,width=6)
# pal = brewer.pal(length($seq2),'Spectral')
# plot($(degsort[1]),fitted(M_t[[1]]),xlim=c(1,10),ylim=c(1,8),col=pal[1],type='l',lwd=2,xlab='Degree',ylab='Trophic level')
# """
# for i=2:length(seq2)
#     R"""
#     lines($(degsort[i]),fitted(M_t[[$(i)]]),col=pal[$(i)],lwd=2)
#     sampsize = min(length($(degsort[i])),1000)
#     samp = sample(seq(1:length($(degsort[i]))),size=sampsize)
#     points(jitter($(degsort[i])[samp]),jitter($(trophsort[i])[samp]),pch='.',col=pal[$(i)])
#     """
# end
# R"""
# lines($Pdegsort,fitted(mpool),lwd=2)
# samp = sample(seq(1:length($Pdegsort)),size=1000)
# points(jitter($(Pdegsort)[samp]),jitter($(Ptrophsort)[samp]),pch='.')
# legend(x=180,y=10.2,legend = c('Pool',$(seq[seq2])),col=c('black',pal),lty=1,lwd=2,title='Time',cex=0.7,bty='n')
# dev.off()
# """



#Species richness vs. connectance
bins = [5;10;25;50;100;200;1000;2000;4000];
seq2 = indexin(bins,seq);
speciesrichness = Array{Int64}(undef,reps,length(seq2));
connectance = Array{Float64}(undef,reps,length(seq2));
for i=1:length(seq2)
    t = seq2[i];
    speciesrichness[:,i] = sprich[:,t];
    connectance[:,i] = conn[:,t];
end
# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/richconn.pdf");
namespace = string("$(homedir())/2014_Lego/Enigma/figures/yog/richconn.pdf");
R"""
pdf($namespace,height=10,width=15)
par(mfrow=c(2,4))
#y <- log($(speciesrichness[:,1]))
#x <- log($(connectance[:,1]))
#m <- lm(y ~ x)
plot($(connectance[:,1]),$(speciesrichness[:,1]),pch='.',log='xy',ylim=c(10,200),xlim=c(0.001,0.02),main=paste(c('t=',$(bins[1])),sep=''),xlab='Connectance',ylab='Richness')
#abline(m)
"""
for i=2:length(seq2)
    R"""
    #y <- log($(speciesrichness[:,1]))
    #x <- log($(connectance[:,1]))
    #m <- lm(y ~ x)
    plot($(connectance[:,i]),$(speciesrichness[:,i]),pch='.',log='xy',ylim=c(10,200),xlim=c(0.001,0.02),main=paste(c('t=',$(bins[i])),sep=''),xlab='Connectance',ylab='Richness')
    #abline(m)
    """
end
R"dev.off()"

# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/richconn2.pdf");
namespace = string("$(homedir())/2014_Lego/Enigma/figures/yog/richconn2.pdf");
R"""
pdf($namespace,height=6,width=12)
par(mfrow=c(1,2))
plot($(connectance[:,5]),$(speciesrichness[:,5]),pch='.',ylim=c(50,200),xlim=c(0.002,0.006),main=paste(c('t=',$(bins[5])),sep=''),xlab='Connectance',ylab='Richness')
plot($(connectance[:,9]),$(speciesrichness[:,9]),pch='.',ylim=c(50,200),xlim=c(0.002,0.006),main=paste(c('t=',$(bins[9])),sep=''),xlab='Connectance',ylab='Richness')
#abline(m)
dev.off()
"""





#Extinction and colonization rates
#Number of engineers vs. size of extinction cascade

# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/sim_settings.jld");
d1 = load(namespace);
reps = d1["reps"];
S = d1["S"];
maxits = d1["maxits"];

lcdf = 500;
EXTCDF = SharedArray{Int64}(reps,lcdf);
extratevec = SharedArray{Float64}(reps,lcdf);
engineers = SharedArray{Int64}(reps,maxits);
sprich = SharedArray{Int64}(reps,maxits);
rich = SharedArray{Int64}(reps,maxits);
@sync @parallel for r=1:reps
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    # namespace_rep = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    namespace_rep = string("$(homedir())/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    # namespace_cid = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    namespace_cid = string("$(homedir())/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
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
        sprich[r,t] = sum(CID[1:S,t]);
        rich[r,t] = sum(CID[:,t]);
        #how many engineers?
        spcid = findall(isodd,CID[1:S,t]);
        engineers[r,t] = sum(sum(m_b[spcid,:],2) .> 0);
    end
    
    spdiff = diff(sprich[r,:]);
    rdiff = diff(rich[r,:]);
    
    colpos = findall(x->x>0,spdiff);
    extpos = findall(x->x<0,spdiff);
    extinctions = spdiff[extpos].*-1;
    colonizations = spdiff[colpos];
    
    extrate = extinctions ./ dt[extpos];
    colrate = colonizations ./ dt[colpos];
    
    # extratevec = collect(0:0.0001:maximum(extrate));
    extratevec[r,:] = collect(range(0,maximum(extrate)/lcdf,lcdf));
    
    extcdf = Array{Int64}(undef,lcdf);
    for i=1:lcdf
        extcdf[i] = length(findall(x->x<extratevec[r,i],extrate))
    end
    
    EXTCDF[r,:] = extcdf;    
    
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
EXTCDFpr = Array{Float64}(undef,reps,lcdf);
for i=1:reps
    EXTCDFpr[i,:] = Array(EXTCDF[i,:])/maximum(Array(EXTCDF[i,:]));
end


sortalg = engsort;
# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/engcdf.pdf");
namespace = string("$(homedir())/2014_Lego/Enigma/figures/yog/engcdf.pdf");
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
# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/objcdf.pdf");
namespace = string("$(homedir())/2014_Lego/Enigma/figures/yog/objcdf.pdf");
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


#POTENTIAL COLONIZERS
sprich = SharedArray{Int64}(reps,maxits);
rich = SharedArray{Int64}(reps,maxits);
pc = SharedArray{Int64}(reps,maxits);
@sync @parallel for r=1:reps
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    # namespace_rep = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    namespace_rep = string("$(homedir())/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    # namespace_cid = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    namespace_cid = string("$(homedir())/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    d3 = load(namespace_cid);
    CID = d3["CID"];
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    #Analysis
    for t = 1:maxits
        cid = findall(isodd,CID[:,t]);
        sprich[r,t] = sum(CID[1:S,t]);
        rich[r,t] = sum(CID[:,t]);
        pc[r,t] = length(potcol(sp_v,int_id,cid,a_b,n_b0,athresh,nthresh));   
    end
end
namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/potcol.jld");
save(namespace,
"sprich", sprich,
"rich", rich,
"pc", pc
);

mpc = vec(mean(pc,1));
sdpc = vec(std(pc,1));
propss = vec(mean(sprich,1)) ./ mean(sprich[:,maxits-100:maxits]);

# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/potcol.pdf");
namespace = string("$(homedir())/2014_Lego/Enigma/figures/yog/potcol.pdf");
R"""
library(RColorBrewer)
pdf($namespace,width=8,height=6)
pal = brewer.pal(3,'Set1');
plot($propss,$mpc/$S,type='l',col=pal[1],lwd=3,xlab='Proportion filled',ylab='Available niche space',ylim=c(0.15,0.25))
polygon(x=c($propss,rev($propss)),y=c(($mpc-$sdpc)/$S,rev(($mpc+$sdpc)/$S)),col=paste(pal[1],50,sep=''),border=NA)
# polygon(x=c(seq(1,$maxits),rev(seq(1,$maxits))),y=c(($mpc-$sdpc)/$S,rev(($mpc+$sdpc)/$S)),col=paste(pal[1],50,sep=''),border=NA)
lines($propss,$mpc/$S,col=pal[1],lwd=3)
dev.off()
"""
