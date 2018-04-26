# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");
loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");

namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/sim_settings.jld");
d1 = load(namespace);
reps = d1["reps"];
S = d1["S"];
tmax = d1["tmax"];
# a_thresh = d1["a_thresh"];
# n_thresh = d1["n_thresh"];
# extinctions = d1["extinctions"];
# probs = d1["probs"]; 


# #We need to 'grab' individual cid_r's within the parallel loop
# #Export them here and then import them separately under parallel rep loop
# @time for r=1:reps
#     namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
#     # namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
#     CID = cid_r[r,:,:];
#     save(namespace,
#     "CID", CID);
# end

#define cid_r everywhere


# NOTE: just analyze the timesteps necessary for the sequence

#We want to catch the full assembly early on, and then skip ahead
seq = [collect(2:50);100;200;1000;2000;];
tseqmax = length(seq);

rich = SharedArray{Int64}(reps,tseqmax);
sprich = SharedArray{Int64}(reps,tseqmax);
turnover = SharedArray{Float64}(reps,tseqmax);
res_overlap = SharedArray{Float64}(reps,tseqmax);
conn = SharedArray{Float64}(reps,tseqmax);
conn_ind = SharedArray{Float64}(reps,tseqmax);
avgdegree = SharedArray{Float64}(reps,tseqmax);
res_overlap_dist = SharedArray{Float64}(reps,tseqmax,S);
degrees = SharedArray{Int64}(reps,tseqmax,S);
trophic = SharedArray{Float64}(reps,tseqmax,S);

@sync @parallel for r=1:reps
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    namespace_rep = string("$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    namespace_cid = string("$(homedir())/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
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
    for t = 1:tseqmax
        
        #construct
        tstep = seq[t];
        cid = find(isodd,CID[tstep,:]);
        cid_old = find(isodd,CID[tstep-1,:]); #because we have this, seq can't start at t=1;
        
        rich[r,t], sprich[r,t], turnover[r,t], res_overlap[r,t], res_overlap_all, conn[r,t], conn_ind[r,t] = dynstructure(cid,cid_old,sp_v,a_b,tp_m,tind_m);     
        
        res_overlap_dist[r,t,1:length(res_overlap_all)] = res_overlap_all; 
        
        deg,troph = structure(S,cid,sp_v,tind_m);
        
        degrees[r,t,1:length(deg)] = deg;
        trophic[r,t,1:length(troph)] = troph;
        avgdegree[r,t] = mean(degrees[r,t,1:length(deg)]);
        
        
    end

end



h = load(string("$(homedir())/2014_Lego/Anime/data/intm_structure.jld"));
Pconn = h["Pconn"];
Pconn_ind = h["Pconn_ind"];
Pres_overlap_dist = h["Pres_overlap_dist"];
Pdegrees = h["Pdegrees"];
Ptl = h["Ptl"];
Preps = size(Pconn)[1];





############
#Connectance
############

lfseq = find(x->x>1,diff(seq))[1];
#This is the initial assembly process
init_conn = conn[:,1:lfseq];
init_conn_ind = conn_ind[:,1:lfseq];
init_conn_trim = Array{Float64}(reps,lfseq)*0;
init_conn_ind_trim = Array{Float64}(reps,lfseq)*0;
for r=1:reps
    init_conn_rm = init_conn[r,find(!iszero,init_conn[r,:])];
    init_conn_ind_rm = init_conn_ind[r,find(!iszero,init_conn_ind[r,:])];
    init_conn_trim[r,1:length(init_conn_rm)] = init_conn_rm;
    init_conn_ind_trim[r,1:length(init_conn_ind_rm)] = init_conn_ind_rm;
end
bins = [2;5;10;50;100;1000;2000];
initsteps = bins[bins.<lfseq]; #use these locations for init
laststeps = bins[bins.>=lfseq]; #use these locations for the rest
lastbins = indexin(laststeps,seq);

#Stitch together
seq_stitch = [initsteps;laststeps];
conn_stitch = [init_conn_trim[:,initsteps] conn[:,lastbins]];

namespace = string("$(homedir())/2014_Lego/Anime/figures2/conn_time2.pdf");
R"""
pdf($namespace,height=5,width=6)
boxplot($(conn_stitch),ylim=c(0,0.1),outline=FALSE,names=$(seq_stitch),
xlab='Time',ylab='Connectance',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(vec(mapslices(mean,conn_stitch,1))),ylim=c(0,0.1),pch=16)
lines($(vec(mapslices(mean,conn_stitch,1))),ylim=c(0,0.1),lwd=2)
lines(seq(0.001,3000),rep(mean($Pconn),3000),lty=2)
dev.off()
"""


############
#Resource overlap
############
#Species pool
Preps = size(Pres_overlap_dist)[1];
Pmeanoverlap = Array{Float64}(Preps)
for r=1:Preps
    Pmeanoverlap[r] = mean(Pres_overlap_dist[r,!isnan(Pres_overlap_dist[r,:])]);
end

lfseq = find(x->x>1,diff(seq))[1];
#This is the initial assembly process
init_overlap = res_overlap[:,1:lfseq];
init_overlap_trim = Array{Float64}(reps,lfseq)*0;
for r=1:reps
    init_overlap_rm = init_overlap[r,find(x->x>0,init_overlap[r,:])];
    init_overlap_trim[r,1:length(init_overlap_rm)] = init_overlap_rm;
end
bins = [2;5;10;50;100;1000;2000];
initsteps = bins[bins.<lfseq]; #use these locations for init
laststeps = bins[bins.>=lfseq]; #use these locations for the rest
lastbins = indexin(laststeps,seq);

#Stitch together
seq_stitch = [initsteps;laststeps];
overlap_stitch = [init_overlap_trim[:,initsteps] res_overlap[:,lastbins]];

namespace = string("$(homedir())/2014_Lego/Anime/figures2/trophicoverlap_time.pdf");
R"""
pdf($namespace,height=5,width=6)
boxplot($(overlap_stitch),ylim=c(0,0.1),outline=FALSE,names=$(seq_stitch),
xlab='Time',ylab='Resource overlap',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(vec(mapslices(mean,overlap_stitch,1))),ylim=c(0,0.1),pch=16)
lines($(vec(mapslices(mean,overlap_stitch,1))),ylim=c(0,0.1),lwd=2)
lines(seq(0.001,3000),rep(mean($Pmeanoverlap),3000),lty=2)
dev.off()
"""

#######################
# DEGREE DISTRIBUTION
#######################


# Somehow display how the *realized* degree distribution changes
bins = [5;10;25;50;100;200;1000;2000;];
seq2 = indexin(bins,seq);

mdegt = Array{Float64}(length(seq2),S)*0;
sddegt = Array{Float64}(length(seq2),S)*0;
t_tic = 0;
for t = seq2
    t_tic = t_tic + 1;
    deg = degrees[:,t,:];
    #Sort each row
    degsort = sort(deg,2,rev=true);
    #which column is the last non-zero?
    lastcol = find(iszero,sum(degsort,1))[1]-1;
    mdeg = Array{Float64}(lastcol);
    sddeg = Array{Float64}(lastcol);
    #Take means but ignore zeros for each column through lascol
    for i=1:lastcol
        mdeg[i] = mean(degsort[!iszero.(degsort[:,i]),i]);
        sddeg[i] = std(degsort[!iszero.(degsort[:,i]),i]);
    end
    mdegt[t_tic,1:length(mdeg)]=mdeg;
    sddegt[t_tic,1:length(mdeg)]=sddeg;
end
Pdegreesort = sort(Pdegrees,2,rev=true);
Pmeandegree = vec(mapslices(mean,Pdegreesort,1));
Psddeg = vec(mapslices(std,Pdegreesort,1));

meanrich = convert(Int64,round(mean(sprich[:,tseqmax]),0));
Pdegreesort = Array{Int64}(5000,meanrich)*0;
for i=1:Preps
    Pdegreesort[i,:] = sort(sample(Pdegrees[i,:],meanrich),rev=true);
end
Pmeandegree = vec(mapslices(mean,Pdegreesort,1));
firstone = find(x->x==1,Pmeandegree)[1];
Pmeandegree[firstone:length(Pmeandegree)] = 1;
Psddeg = vec(mapslices(std,Pdegreesort,1));
Psddeg[firstone:length(Psddeg)] = 0;


namespace = string("$(homedir())/2014_Lego/Anime/figures2/degreedist_time2.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
pal = brewer.pal($(length(seq2)),'Spectral')
numsp = length($(mdegt[1,!iszero.(mdegt[1,:])]))
plot($(mdegt[1,!iszero.(mdegt[1,:])]),xlim=c(1,250),ylim=c(1,100),log='y',col=pal[1],type='l',lwd=2,xlab = 'Number of species', ylab='Mean degree')
sdev_pre = $(sddegt[1,find(x->x>0,sddegt[1,:])]);
sdev = numeric(length($(mdegt[1,find(x->x>0,mdegt[1,:])])))
sdev[1:length(sdev_pre)]=sdev_pre
polygon(x=c(seq(1,length(sdev)),seq(length(sdev),1)),
y=c($(mdegt[1,!iszero.(mdegt[1,:])])[1:length(sdev)]+sdev,
rev($(mdegt[1,!iszero.(mdegt[1,:])])[1:length(sdev)]-sdev)),col=paste(pal[1],65,sep=''),border=NA)
lines($(mdegt[1,!iszero.(mdegt[1,:])]),xlim=c(1,200),ylim=c(0.01,50),log='y',col=pal[1],type='l',lwd=2,xlab = 'Number of species', ylab='Median degree')
"""
for i=2:length(seq2)
    R"""
    sdev_pre = $(sddegt[i,find(x->x>0,sddegt[i,:])]);
    sdev = numeric(length($(mdegt[i,find(x->x>0,mdegt[i,:])])))
    sdev[1:length(sdev_pre)]=sdev_pre
    polygon(x=c(seq(1,length(sdev)),seq(length(sdev),1)),
    y=c($(mdegt[i,!iszero.(mdegt[i,:])])[1:length(sdev)]+sdev,
    rev($(mdegt[i,!iszero.(mdegt[i,:])])[1:length(sdev)]-sdev)),col=paste(pal[$i],65,sep=''),border=NA)
    lines($(mdegt[i,!iszero.(mdegt[i,:])]),col=pal[$i],lwd=2)
    """
end
R"""
maxsp = $meanrich;
sdev = $(Psddeg)[1:maxsp];
polygon(x=c(seq(1,maxsp),seq(maxsp,1)),
y=c($(Pmeandegree)[1:maxsp]+sdev[1:maxsp],
rev($(Pmeandegree)[1:maxsp]-sdev[1:maxsp])),col=paste('#000000',65,sep=''),border=NA)
lines($(Pmeandegree)[1:maxsp],col='black',type='l',lwd=2)

legend(x=215,y=120,legend = c('Pool',$(seq[seq2])),col=c('black',pal),lty=1,lwd=2,title='Time',cex=0.7,bty='n')
dev.off()
"""


#####################
# Trophic levels NOTE: TODO
#####################
bins = [10;25;50;100;200;1000;2000;];
seq2 = indexin(bins,seq);

R"M_t=list()"
degsort = Array{Array{Int64}}(length(seq2));
trophsort = Array{Array{Float64}}(length(seq2));
for i=1:length(seq2)
    t = seq2[i];
    deg = reshape(degrees[:,t,:],reps*S);
    troph = reshape(trophic[:,t,:],reps*S);
    deg_trim = deg[find(!iszero,deg)];
    trophic_trim = troph[find(!iszero,troph)];
    x = deg_trim;
    y = trophic_trim;
    z = [x y];
    sp = sortperm(z[:,1]);
    zsort = [z[sp,1] z[sp,2]];
    xsort = zsort[:,1];
    ysort = zsort[:,2];
    degsort[i] = xsort;
    trophsort[i] = ysort;
    R"""
    y <- $ysort
    x <- $xsort
    cf <- c(0,0,0)
    m <- try(nls(y ~ a + b * I(x^z), start = list(a = 1, b = 1, z = -1)),silent=TRUE)
    M_t[[$i]] = m
    """
end

Pdeg = reshape(Array(Pdegrees),size(Pdegrees)[1]*size(Pdegrees)[2]);
Ptroph = reshape(Array(Ptl),size(Ptl)[1]*size(Ptl)[2]);
attached = find(!iszero,Pdeg);
Pdeg = Pdeg[attached];
Ptroph = Ptroph[attached];
z = [Pdeg Ptroph];
sp = sortperm(z[:,1]);
zsort = [z[sp,1] z[sp,2]];
xsort = zsort[:,1];
ysort = zsort[:,2];
Pdegsort = xsort;
Ptrophsort = ysort;
R"""
y <- $ysort
x <- $xsort
mpool <- try(nls(y ~ a + b * I(x^z), start = list(a = 1, b = 1, z = -1)),silent=TRUE)
"""


namespace = string("$(homedir())/2014_Lego/Anime/figures2/trophic_degrees_time.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
pal = brewer.pal(length($seq2),'Spectral')
plot($(degsort[1]),fitted(M_t[[1]]),xlim=c(1,400),ylim=c(1,10),log='x',col=pal[1],type='l',lwd=2,xlab='Degree',ylab='Trophic level')
"""
for i=2:length(seq2)
    R"""
    lines($(degsort[i]),fitted(M_t[[$i]]),col=pal[$i],lwd=2)
    sampsize = min(length($(degsort[i])),1000)
    samp = sample(seq(1:length($(degsort[i]))),size=sampsize)
    points(jitter($(degsort[i])[samp]),jitter($(trophsort[i])[samp]),pch='.',col=pal[$i])
    """
end
R"""
lines($Pdegsort,fitted(mpool),lwd=2)
samp = sample(seq(1:length($Pdegsort)),size=1000)
points(jitter($(Pdegsort)[samp]),jitter($(Ptrophsort)[samp]),pch='.')
legend(x=180,y=10.2,legend = c('Pool',$(seq[seq2])),col=c('black',pal),lty=1,lwd=2,title='Time',cex=0.7,bty='n')
dev.off()
"""

