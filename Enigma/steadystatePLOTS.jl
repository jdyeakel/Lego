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
sprichinweb = SharedArray{Int64}(reps,tseqmax);
sprichinwebnoprim = SharedArray{Int64}(reps,tseqmax);

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

realizeddegree = SharedArray{Float64}(reps,tseqmax,S);
potentialdegree = SharedArray{Float64}(reps,tseqmax,S);

realizedG = SharedArray{Float64}(reps,tseqmax,S);
realizedGnoprim = SharedArray{Float64}(reps,tseqmax,S);
realizedGavgc = SharedArray{Float64}(reps,tseqmax,S);

potentialG = SharedArray{Float64}(reps,tseqmax,S);
C = SharedArray{Float64}(reps,tseqmax);

#Save richness for colonizing phase
coltraj = SharedArray{Int64}(reps,1000);
clocktraj = SharedArray{Float64}(reps,1000);

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
    
    #Save FULL colonizing phase
    coltraj[r,:] = sum(CID[:,1:1000],dims=1);
    clocktraj[r,:] = clock[1:1000];

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
        
        # realizeddegree[r,t] = mean(vec(sum(a_b[cid,[1;cid]],dims=2)));
        # potentialdegree[r,t] = mean(vec(sum(a_b[cid,:],dims=2)));
        realizeddegree[r,t,1:length(cid)] = vec(sum(a_b[cid,[1;cid]],dims=2));
        potentialdegree[r,t,1:length(cid)] = vec(sum(a_b[cid,:],dims=2));
        
        
        #NOTE species disconnected from web can be counted
        
        #this are species that are connected to the foodweb in any way
        cidinweb = sort(cid[findall(!iszero,vec(sum(a_b[cid,[1;cid]],dims=2)))]);
        
        sprichinweb[r,t] = length(cidinweb);
        
        C[r,t] = (sum(a_b[cidinweb,[1;cid]]))/((length(cidinweb))^2);
        
        realizedG[r,t,1:length(cidinweb)] = vec(sum(a_b[cidinweb,[1;cid]],dims=2)) .* (1/((sum(a_b[cidinweb,[1;cid]]))/(length(cidinweb))));
        
        potentialG[r,t,1:length(cidinweb)] = vec(sum(a_b[cidinweb,:],dims=2)) .* (1/((sum(a_b[cidinweb,[1;cid]]))/(length(cidinweb))));
        
        cidinwebnoprim = sort(cid[findall(!iszero,vec(sum(a_b[cid,cid],dims=2)))]);
        sprichinwebnoprim[r,t] = length(cidinwebnoprim);
        
        if sprichinwebnoprim[r,t] > 0
            realizedGnoprim[r,t,1:length(cidinwebnoprim)] = vec(sum(a_b[cidinwebnoprim,cid],dims=2)) .* (1/((sum(a_b[cidinwebnoprim,cid]))/(length(cidinwebnoprim))));
        end
        
        #compared to steady state C averaged across webs
        #AS CLOSE TO PIETCHNICK AS POSSIBLE
        SSC = 0.0108;
        if sprichinwebnoprim[r,t] > 0
            realizedGavgc[r,t,1:length(cidinweb)] = vec(sum(a_b[cidinweb,[1;cid]],dims=2)) .* (1/(SSC*length(cidinweb)));
        end
        
    end

end

filename = "figures/logistic_test.pdf";
namespace = smartpath(filename);
#logistic fit
R"""
d <- data.frame("Time"=$(clocktraj[1,:]),"Rich"=$(coltraj[1,:]))

Time = as.numeric(d["Time"][,1]);
Rich = as.numeric(d["Rich"][,1]);
SS<-getInitial(Rich~SSlogis(Time,alpha,xmid,scale),data=d)
K_start<-SS["alpha"]
R_start<-1/SS["scale"]
N0_start<-SS["alpha"]/(exp(SS["xmid"]/SS["scale"])+1)
log_formula<-formula(Rich~K*N0*exp(R*Time)/(K+N0*(exp(R*Time)-1)))
m<-nls(log_formula,start=list(K=K_start,R=R_start,N0=N0_start))

pdf($namespace,width=6,height=5)
plot($(clocktraj[1,:]),$(coltraj[1,:]))
curve(predict(m,data.frame(Time=x),type="resp"),add=TRUE)
dev.off()
"""



filename = "data/intm_structure.jld";
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
bins = [2;5;10;20;30;40;100;1000;4000];
conn_stitch,seq_stitch = sortassembly(conn,bins,seq);
meanconn = [mean(conn_stitch[findall(!isnan,conn_stitch[:,i]),i]) for i=1:length(seq_stitch)];
sreps = sample(collect(1:reps),100);



rsample=sample(collect(1:reps),100);
r=rsample[1];
filename = "data/steadystate/int_m.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace int_m tp_m tind_m mp_m mind_m;
filename = "data/steadystate/cid.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace CID clock;
filename = "figures/yog/conn_time2.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1');
pdf($namespace,height=5,width=10)
par(mfrow=c(1,2))
plot($(collect(1:4000)),$(vec(sum(CID,dims=1))),type='l',log='x',col=paste(pal[2],40,sep=''),xlim=c(1,4000),ylim=c(0,180),xlab='Time',ylab='Species richness')
"""
for i=2:length(rsample)
    r=rsample[i];
    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;
    # filename = "figures/yog/assembly_time.pdf"
    # namespace = smartpath(filename);
    R"""
    lines($(collect(1:4000)),$(vec(sum(CID,dims=1))),type='l',col=paste(pal[2],40,sep=''))
    """
end
R"""
plot(jitter(rep(1,$reps)),jitter($(conn_stitch[:,1])),pch='.',col=pal[2],xlim=c(1,length($seq_stitch)),ylim=c(0.0,0.2),xaxt='n',
xlab='Time',ylab=expression(Connectance ~~ L/S^{2}))
"""
for t=2:length(seq_stitch)
    R"""
    points(jitter(rep($(t),$reps)),jitter($(conn_stitch[:,t])),pch='.',col=pal[2])
    """
end
R"""
boxplot($(conn_stitch),outline=FALSE,names=$(seq_stitch),xlab='',ylab='',
pars = list(boxwex = 0.2, staplewex = 0.5, outwex = 0.5,whisklty=1),add=TRUE)
points($(meanconn),ylim=c(0,0.03),pch=16)
lines($(meanconn),ylim=c(0,0.03),lwd=2)
lines(seq(0.001,3000),rep(mean($Pconn),3000),lty=3)
dev.off()
"""



filename = "figures/yog/conn_time.pdf"
namespace = smartpath(filename);

R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,height=5,width=10)
par(mfrow=c(1,2))
plot($(sprich[sreps])
plot(jitter(rep(1,$reps)),jitter($(conn_stitch[:,1])),pch='.',col=pal[2],xlim=c(1,length($seq_stitch)),ylim=c(0.0,0.2),xaxt='n',
xlab='Time',ylab=expression(Connectance ~~ L/S^{2}))
"""
for t=2:length(seq_stitch)
    R"""
    points(jitter(rep($(t),$reps)),jitter($(conn_stitch[:,t])),pch='.',col=pal[2])
    """
end
R"""
boxplot($(conn_stitch),outline=FALSE,names=$(seq_stitch),xlab='',ylab='',
pars = list(boxwex = 0.2, staplewex = 0.5, outwex = 0.5,whisklty=1),add=TRUE)
points($(meanconn),ylim=c(0,0.03),pch=16)
lines($(meanconn),ylim=c(0,0.03),lwd=2)
lines(seq(0.001,3000),rep(mean($Pconn),3000),lty=3)
dev.off()
"""


#mutualistic connectance
mutconn_stitch,seq_stitch = sortassembly(mutconn,bins,seq);
meanmutconn = [mean(mutconn_stitch[findall(!isnan,mutconn_stitch[:,i]),i]) for i=1:length(seq_stitch)]


filename = "figures/yog/mutconn_time.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=6)
boxplot($(mutconn_stitch),ylim=c(0.0001,0.01),log='y',outline=FALSE,names=$(seq_stitch),
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
meanoverlap = [mean(overlap_stitch[isnan.(overlap_stitch[:,i]).==false,i]) for i=1:length(seq_stitch)];

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
meanoverlap = [mean(useroverlap_stitch[isnan.(useroverlap_stitch[:,i]).==false,i]) for i=1:length(seq_stitch)];
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

filename = "figures/yog/degreedist_time2.pdf";
namespace = smartpath(filename);
# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/yog/degreedist_time2.pdf");
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


#Real degreea and trophic distribution

##################
#USE THESE 6/24/2019
##################

bins = [5;10;25;50;100;200;500;1000;2000;4000;];
# bins = [5;50;100;4000;];
seq2 = indexin(bins,seq);
tmaxdegree = zeros(Int64,length(seq2));
tmaxtrophic = zeros(Int64,length(seq2));
maxdegree = Array{Int64}(undef,reps,length(seq2));
maxtrophic = Array{Int64}(undef,reps,length(seq2));
freqdegreereps = Array{Array}(undef,reps);
freqtrophicreps = Array{Array}(undef,reps);
for r=1:reps
    freqdegreetime = Array{Array}(undef,length(seq2));
    freqtrophictime = Array{Array}(undef,length(seq2));
    let tic = 0
        for t=seq2
            tic += 1;
            alldegrees = degrees[r,t,:][degrees[r,t,:] .> 0];
            alltrophic = trophic[r,t,:][trophic[r,t,:] .>= 0.9];
            maxdegree[r,tic] = maximum(alldegrees);
            maxtrophic[r,tic] = round(Int64,maximum(alltrophic))+1;
            freqdegree = Array{Float64}(undef,maxdegree[r,tic]);
            freqtrophic = Array{Float64}(undef,maxtrophic[r,tic]+1);
            for i=1:maxdegree[r,tic]
                freqdegree[i] = length(findall(x->x==i,alldegrees))/length(alldegrees);
            end
            for i=0:maxtrophic[r,tic]
                freqtrophic[i+1] = length(findall(x->(x>=i && x<i+1),alltrophic))/length(alltrophic);
            end
            freqdegreetime[tic] = freqdegree;
            freqtrophictime[tic] = freqtrophic;
            tmaxdegree[tic] = maximum([maxdegree[r,tic],tmaxdegree[tic]])
            tmaxtrophic[tic] = maximum([maxtrophic[r,tic],tmaxtrophic[tic]])
        end
    end
    freqdegreereps[r] = freqdegreetime;
    freqtrophicreps[r] = freqtrophictime;
end
#lets ignore results with >20 trophic levels!
toignore = findall(!iszero,vec(sum(maxtrophic .> 20,dims=2)));
newreps = setdiff(collect(1:reps),toignore);


maxdegreeall = maximum(tmaxdegree);
# maxtrophic = maximum(tmaxtrophic)+1;
maxtrophicall = 21;
#Reorganizing
degreedistreps = zeros(Float64,length(newreps),length(seq2),maxdegreeall);
trophicdistreps = zeros(Float64,length(newreps),length(seq2),maxtrophicall);
let rtic=0
    for r=newreps
        rtic += 1;
        for t=1:length(seq2)
            degreedistreps[rtic,t,1:length(freqdegreereps[r][t])] = freqdegreereps[r][t];
            trophicdistreps[rtic,t,1:length(freqtrophicreps[r][t])] = freqtrophicreps[r][t];
        end
    end
end
#means over reps
meandegreedist = Array{Float64}(undef,length(seq2),maxdegreeall);
meantrophicdist = Array{Float64}(undef,length(seq2),maxtrophicall);
for t=1:length(seq2)
    meandegreedist[t,:] = mean(degreedistreps[:,t,:],dims=1);
    meantrophicdist[t,:] = mean(trophicdistreps[:,t,:],dims=1);
end


filename = "figures/yog/degreedist_time3trimunlogged.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($seq2))
pdf($namespace,height=5,width=6)
plot(seq(1,length($(meandegreedist[1,:]))),$(meandegreedist[1,:]),type='l',xlim=c(1,10),ylim=c(0.000001,1),col=pal[1],xlab='Species degree',ylab='Probability')
points(seq(1,length($(meandegreedist[1,:]))),$(meandegreedist[1,:]),pch=21,bg=pal[1],col='black')
"""
for i=2:length(seq2)
    R"""
    lines(seq(1,length($(meandegreedist[i,:]))),$(meandegreedist[i,:]),col=pal[$i])
    points(seq(1,length($(meandegreedist[i,:]))),$(meandegreedist[i,:]),pch=21,bg=pal[$i],col='black')
    """
end
R"""
legend(x=8,y=1,legend=$(seq[seq2]),col=pal,cex=0.8,pch=16,bty='n',title='Assembly time')
dev.off()
"""




filename = "figures/yog/degreedist_time3trim.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($seq2))
pdf($namespace,height=5,width=6)
plot(seq(1,length($(meandegreedist[1,:]))),$(meandegreedist[1,:]),type='l',xlim=c(1,10),ylim=c(0.000001,1),col=pal[1],xlab='Species degree',ylab='Probability',log='y')
points(seq(1,length($(meandegreedist[1,:]))),$(meandegreedist[1,:]),pch=21,bg=pal[1],col='black')
"""
for i=2:length(seq2)
    R"""
    lines(seq(1,length($(meandegreedist[i,:]))),$(meandegreedist[i,:]),col=pal[$i])
    points(seq(1,length($(meandegreedist[i,:]))),$(meandegreedist[i,:]),pch=21,bg=pal[$i],col='black')
    """
end
R"""
legend(x=8,y=1,legend=$(seq[seq2]),col=pal,cex=0.8,pch=16,bty='n',title='Assembly time')
dev.off()
"""



filename = "figures/yog/trophicdist_time3.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal($(length(seq2)),'Spectral')
pdf($namespace,height=5,width=6)
plot(seq(1,length($(meantrophicdist[1,:]))),$(meantrophicdist[1,:]),type='l',xlim=c(1,10),ylim=c(0.000000,0.8),col=pal[1],xlab='Trophic level',ylab='Probability')
points(seq(1,length($(meantrophicdist[1,:]))),$(meantrophicdist[1,:]),pch=21,bg=pal[1],col='black')
"""
for i=2:length(seq2)
    R"""
    lines(seq(1,length($(meantrophicdist[i,:]))),$(meantrophicdist[i,:]),col=pal[$i])
    points(seq(1,length($(meantrophicdist[i,:]))),$(meantrophicdist[i,:]),pch=21,bg=pal[$i],col='black')
    """
end
R"""
legend(x=8,y=0.8,legend=$(seq[seq2]),col=colorRampPalette(brewer.pal(9,"Spectral"))(9),cex=0.8,pch=16,bty='n',title='Assembly time')
dev.off()
"""



filename = "../manuscript/fig_trophic.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($seq2))
pdf($namespace,height=5,width=6)
fulldist = $(meantrophicdist[1,:]);
trimdist = fulldist[which(fulldist>0.005)];
plot(trimdist,seq(1,length(trimdist)),type='l',xlim=c(0,0.6),ylim=c(1,12),col=pal[1],xlab='Frequency',ylab='Trophic level')
points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[1],col='black')
"""
for i=1:12
    R"""
    rect(0,$i-0.5,$(meantrophicdist[length(seq2),i]),$i+0.5,col=paste(pal[length($seq2)],50,sep=''),border=NA)
    """
end
for i=1:length(seq2)
    R"""
    fulldist = $(meantrophicdist[i,:]);
    trimdist = fulldist[which(fulldist>0.005)];
    lines(trimdist,seq(1,length(trimdist)),col=pal[$i])
    points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[$i],col='black')
    """
end
R"""
legend(x=0.45,y=12,legend=$(seq[seq2]),col=pal,cex=0.8,pch=16,bty='n',title='Assembly time')
dev.off()
"""



#Resource specialization of realized vs. potential niche over time

bins = [5;10;25;50;100;200;500;1000;2000;4000;];
# bins=seq;
seq2 = indexin(bins,seq);
meanpotentialdegree = Array{Float64}(undef,reps,length(seq2));
meanrealizeddegree = Array{Float64}(undef,reps,length(seq2));
rG = Array{Float64}(undef,reps,length(seq2));
pG = Array{Float64}(undef,reps,length(seq2));
propG = Array{Float64}(undef,reps,length(seq2));
propGavgc = Array{Float64}(undef,reps,length(seq2));
propGnoprim = Array{Float64}(undef,reps,length(seq2));
uppertrophicpropG = Array{Float64}(undef,reps,length(seq2));
potpropG = Array{Float64}(undef,reps,length(seq2));
meantrophic = Array{Float64}(undef,reps,length(seq2));
meandegrees = Array{Float64}(undef,reps,length(seq2));

for r=1:reps
    for t=1:length(seq2)
        
        meanpotentialdegree[r,t] = mean(potentialdegree[r,seq2[t],:][potentialdegree[r,seq2[t],:].>0]);
        meanrealizeddegree[r,t] = mean(realizeddegree[r,seq2[t],:][realizeddegree[r,seq2[t],:].>0]);
        
        rG[r,t] = mean(realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0]);
        pG[r,t] = mean(potentialG[r,seq2[t],:][potentialG[r,seq2[t],:].>0]);
        
        # propG[r,t] = sum((realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0] .* C[r,seq2[t]]) .> (C[r,seq2[t]]))/sum(realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0]);
        
        # propG[r,t] = sum((realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0] ) .> (1))/sum(realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0]);
        
        propG[r,t] = sum((realizedG[r,seq2[t],1:sprichinweb[r,seq2[t]]]) .> (1))/sprichinweb[r,seq2[t]];
        
        propGnoprim[r,t] = sum((realizedGnoprim[r,seq2[t],1:sprichinwebnoprim[r,seq2[t]]]) .> (1))/sprichinwebnoprim[r,seq2[t]];
        
        propGavgc[r,t] = sum((realizedGavgc[r,seq2[t],1:sprichinwebnoprim[r,seq2[t]]]) .> (1))/sprichinwebnoprim[r,seq2[t]];
        
        potpropG[r,t] = sum((potentialG[r,seq2[t],1:sprichinweb[r,seq2[t]]]) .> (1))/sprichinweb[r,seq2[t]];
        
        uppertrophicpropG[r,t] = sum((realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])][realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])].>0]) .> (1))/sum(realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])][realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])].>0]);
        
        meantrophic[r,t] = mean(trophic[r,seq2[t],:][trophic[r,seq2[t],:] .> 0]);
        meandegrees[r,t] = mean(degrees[r,seq2[t],:][degrees[r,seq2[t],:] .> 0]);
        
    end
end

#which reps have low values of propG durin timesteps 1:100?
earlygen = unique([findall(x->x<0.5,1 .- propG[:,3]);findall(x->x<0.5,1 .- propG[:,3])]);
earlysp = setdiff(collect(1:reps),earlygen);

mps = mean(meanpotentialdegree,dims=1);
mrs = mean(meanrealizeddegree,dims=1);


mpG = mean(potpropG,dims=1);
mrG = meanfinite(propG,1);
mrGgen = meanfinite(propG[earlygen,:],1);
mrGnoprim = meanfinite(propGnoprim,1);
mrGavgc = meanfinite(propGavgc,1);
minmrGnoprim = Array{Float64}(undef,length(seq2));
maxmrGnoprim = Array{Float64}(undef,length(seq2));
for i=1:length(seq2)
    minmrGnoprim[i]=quantile(propGnoprim[isnan.(propGnoprim[:,i]).==false,i],0.05);
    maxmrGnoprim[i]=quantile(propGnoprim[isnan.(propGnoprim[:,i]).==false,i],0.95);
end
# mrGuppertrophic = meanfinite(uppertrophicpropG,1);

filename = "figures/yog/specialization.pdf"
namespace = smartpath(filename);
# R"""
# library(RColorBrewer)
# pal = brewer.pal($(length(seq2)),'Spectral')
# pdf($namespace,height=5,width=10)
# par(mfrow=c(1,2))
# plot(jitter($(meanrealizeddegree[:,1])),jitter($(meanpotentialdegree[:,1])),pch='.',col=pal[1],xlim=c(1,2),ylim=c(1,4),xlab='Mean realized trophic links',ylab='Mean potential trophic links')
# # points(jitter($(meanpotentialdegree[earlygen,1])),jitter($(meanrealizeddegree[earlygen,1])),pch='.',col='black')
# """
# for i=2:length(seq2)
#     R"""
#     points(jitter($(meanrealizeddegree[:,i])),jitter($(meanpotentialdegree[:,i])),pch='.',col=pal[$i])
#     # points(jitter($(meanpotentialdegree[earlygen,i])),jitter($(meanrealizeddegree[earlygen,i])),pch='.',col='black')
#     """
# end
# R"""
# legend(x=1.8,y=4.05,legend=$(seq[seq2]),col=colorRampPalette(brewer.pal(9,"Spectral"))(9),cex=0.7,pch=16,bty='n',title='Assembly time')
# lines(seq(-1,5),seq(-1,5))
# lines($mrs,$mps)
# points($mrs,$mps,pch=21,col='black',bg=pal)
R"""
library(RColorBrewer)
pal = brewer.pal($(length(seq2)),'Spectral')
pdf($namespace,height=5,width=6)
par(mfrow=c(1,1))
plot(jitter($(repeat([seq[seq2[1]]],reps))),1-($(propGnoprim[:,1])),pch='.',cex=3,col=pal[1],xlim=c(5,4000),ylim=c(0,1),xlab='Time',ylab='Proportion specialists',log='x',cex.axis=0.85)
"""
for i=2:length(seq2)
    R"""
    points(jitter($(repeat([seq[seq2[i]]],reps))),1-($(propGnoprim[:,i])),pch='.',cex=3,col=pal[$i])
    """
    # nona = 1 .- propGnoprim[isnan.(propGnoprim[:,i]).==false,i];
    # qrs = nona[findall(x->x<quantile(nona,0.05),nona)];
    # R"""
    # points(jitter($(repeat([seq[seq2[i]]],length(qrs)))),$qrs,pch=1,col=pal[$i])
    # """
end
R"""
lines($(seq[seq2]),1-$mrG)
points($(seq[seq2]),1-$mrG,pch=21,col='black',bg=pal,cex=1.5)
# lines($(seq[seq2]),1-$mrGnoprim)
# points($(seq[seq2]),1-$mrGnoprim,pch=22,col='black',bg=pal)
lines($(seq[seq2]),1-$mrGavgc)
points($(seq[seq2]),1-$mrGavgc,pch=23,col='black',bg=pal,cex=1.5)
# lines($(seq[seq2]),1-$mrGgen)
# points($(seq[seq2]),1-$mrGgen,pch=21,col='black',bg=pal)
# lines($(seq[seq2]),1-$mrGuppertrophic)
# points($(seq[seq2]),1-$mrGuppertrophic,pch=23,col='black',bg=pal)
# lines($(seq[seq2]),1-$mpG)
# points($(seq[seq2]),1-$mpG,pch=24,col='black',bg=pal)
# legend(x=1.8,y=4.05,legend=$(seq[seq2]),col=colorRampPalette(brewer.pal(9,"Spectral"))(9),cex=0.7,pch=16,bty='n',title='Assembly time')
dev.off()
"""


#######################
# FOR MANUSCRIPT: top panel specialization; bottom panel trophic
#########################################



filename = "../manuscript/fig_trophic.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pdf($namespace,height=7,width=5)
layout(matrix(c(1,2), 2, 1, byrow = TRUE), 
   widths=c(1,1), heights=c(0.4,0.5))
par(oma = c(0.5, 1, 1, 1), mar = c(3, 4, 0, 1))
#SPECIALIZATION
pal = brewer.pal($(length(seq2)),'Spectral')
timelabels = parse(text=c("5","10","25","50",paste("10","^2"),paste("2.10","^2"),paste("5*10","^2"),paste("10","^3"),paste("2*10","^3"),paste("4*10","^3")))
# par(mfrow=c(2,1))
plot($(seq[seq2]),1 - $mrGavgc, xlim=c(5,4000), ylim=c(0,1), xlab='', ylab='', log='x', cex.axis=0.85, pch=23, col='black', bg=pal, cex=1.5,axes=FALSE)
axis(2,at=seq(0,1,by=0.2),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=$(seq[seq2]),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Proportion specialists', line=2.5, cex.lab=1.2)
title(xlab='Assembly time', line=2.0, cex.lab=1.2)
"""
for i=1:length(seq2)
    R"""
    propspec = 1-($(propGavgc[:,i]))
    propspectrim = propspec[which(propspec > 0)]
    if (length(propspectrim)>0) {
        points(jitter(rep($(seq[seq2[i]]),length(propspectrim))),propspectrim,pch='.',cex=3,col=pal[$i])
    }
    """
end
R"""
lines($(seq[seq2]),1-$mrGavgc,lwd=2)
points($(seq[seq2]),1-$mrGavgc,pch=23,col='black',bg=pal,cex=1.5)

#TROPHIC

pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($seq2))
fulldist = $(meantrophicdist[1,:]);
trimdist = fulldist[which(fulldist>0.005)];
plot(trimdist,seq(1,length(trimdist)),type='l',xlim=c(0,0.6),ylim=c(1,12),col=pal[1],xlab='',ylab='',axes=FALSE)
axis(2,at=seq(1:12),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=seq(0:0.6,by=0.2),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Trophic level (TL)', line=2.5, cex.lab=1.2)
title(xlab='Frequency', line=2.0, cex.lab=1.2)
points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[1],col='black')
"""
for i=1:12
    R"""
    rect(0,$i-0.5,$(meantrophicdist[length(seq2),i]),$i+0.5,col=paste(pal[length($seq2)],50,sep=''),border=NA)
    """
end
for i=1:length(seq2)
    R"""
    fulldist = $(meantrophicdist[i,:]);
    trimdist = fulldist[which(fulldist>0.005)];
    lines(trimdist,seq(1,length(trimdist)),col=pal[$i])
    points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[$i],col='black')
    """
end
R"""
legend(x=0.4,y=12,legend=$(seq[seq2]),col=pal,cex=1,pch=16,bty='n',title='Assembly time')

dev.off()
"""







filename = "../manuscript/fig_trophic2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
   widths=c(1,1,1), heights=c(0.8,1,1))
par(oma = c(0.5, 1, 1, 1), mar = c(3, 4, 1, 1))
"""
rsample=sample(collect(1:reps),100);
r=rsample[1];
filename = "data/steadystate/int_m.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace int_m tp_m tind_m mp_m mind_m;
filename = "data/steadystate/cid.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace CID clock;
R"""
pal=brewer.pal(3,'Set1')
plot($(clock),$(vec(sum(CID,dims=1))),type='l',col=paste(pal[2],40,sep=''),xlim=c(0,100),ylim=c(0,180),xlab='',ylab='',axes=FALSE)
axis(1)
axis(2,las=1)
title(ylab='Species richness', line=3, cex.lab=1.2)
title(xlab='Time', line=2.0, cex.lab=1.2)
mtext(paste0("a"), side = 3, adj = -0.16, 
    line = 0.3,cex=1.2,font=2)
"""

for i=2:length(rsample)
    r=rsample[i];
    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;
    # filename = "figures/yog/assembly_time.pdf"
    # namespace = smartpath(filename);
    R"""
    lines($(clock),$(vec(sum(CID,dims=1))),type='l',col=paste(pal[2],40,sep=''))
    """
end



#SPECIALIZATION
R"""
pal = brewer.pal($(length(seq2)),'Spectral')
timelabels = parse(text=c("5","10","25","50",paste("10","^2"),paste("2.10","^2"),paste("5*10","^2"),paste("10","^3"),paste("2*10","^3"),paste("4*10","^3")))
# par(mfrow=c(2,1))
plot($(seq[seq2]),1 - $mrGavgc, xlim=c(5,4000), ylim=c(0,1), xlab='', ylab='', log='x', cex.axis=0.85, pch=23, col='black', bg=pal, cex=1.5,axes=FALSE)
axis(2,at=seq(0,1,by=0.2),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=$(seq[seq2]),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Proportion specialists', line=3, cex.lab=1.2)
title(xlab='Assembly time', line=2.0, cex.lab=1.2)
mtext(paste0("b"), side = 3, adj = -0.4, 
    line = -0.6,cex=1.2,font=2)
"""
for i=1:length(seq2)
    R"""
    propspec = 1-($(propGavgc[:,i]))
    propspectrim = propspec[which(propspec > 0)]
    if (length(propspectrim)>0) {
        points(jitter(rep($(seq[seq2[i]]),length(propspectrim))),propspectrim,pch='.',cex=3,col=pal[$i])
    }
    """
end
R"""
lines($(seq[seq2]),1-$mrGavgc,lwd=2)
points($(seq[seq2]),1-$mrGavgc,pch=23,col='black',bg=pal,cex=1.5)

#TROPHIC

pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($seq2))
fulldist = $(meantrophicdist[1,:]);
trimdist = fulldist[which(fulldist>0.005)];
plot(trimdist,seq(1,length(trimdist)),type='l',xlim=c(0,0.6),ylim=c(1,12),col=pal[1],xlab='',ylab='',axes=FALSE)
axis(2,at=seq(1:12),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=seq(0:0.6,by=0.2),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Trophic level (TL)', line=2, cex.lab=1.2)
title(xlab='Frequency', line=2.0, cex.lab=1.2)
mtext(paste0("c"), side = 3, adj = -0.3, 
    line = -0.6,cex=1.2,font=2)
points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[1],col='black')
"""
for i=1:12
    R"""
    rect(0,$i-0.5,$(meantrophicdist[length(seq2),i]),$i+0.5,col=paste(pal[length($seq2)],50,sep=''),border=NA)
    """
end
for i=1:length(seq2)
    R"""
    fulldist = $(meantrophicdist[i,:]);
    trimdist = fulldist[which(fulldist>0.005)];
    lines(trimdist,seq(1,length(trimdist)),col=pal[$i])
    points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[$i],col='black')
    """
end
R"""
legend(x=0.45,y=12.5,legend=$(seq[seq2]),pt.cex=1.2,pt.bg=pal,col='black',pch=21,bty='n',title='',cex=0.8)
text(0.45,11.7,'Assembly time')
dev.off()
"""



















filename = "figures/yog/specialization2.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(9,"Spectral"))(length($seq))
pdf($namespace,height=5,width=5)
plot($(seq[seq2]),$mps,type='l',ylim=c(1,3),log='x')
points($(seq[seq2]),$mps,pch=16,col=pal,cex=0.8)
dev.off()
"""


rsample=sample(collect(1:reps),100);
r=rsample[1];
filename = "data/steadystate/int_m.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace int_m tp_m tind_m mp_m mind_m;
filename = "data/steadystate/cid.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace CID clock;
filename = "figures/yog/conn_time2.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1');
pdf($namespace,height=5,width=10)
par(mfrow=c(1,2))
plot($(clock),$(vec(sum(CID,dims=1))),type='l',log='x',col=paste(pal[2],40,sep=''),ylim=c(0,180))
"""
for i=2:length(rsample)
    r=rsample[i];
    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;
    # filename = "figures/yog/assembly_time.pdf"
    # namespace = smartpath(filename);
    R"""
    lines($(clock),$(vec(sum(CID,dims=1))),type='l',col=paste(pal[2],40,sep=''))
    """
end
R"""
plot(jitter(rep(1,$reps)),jitter($(conn_stitch[:,1])),pch='.',col=pal[2],xlim=c(1,length($seq_stitch)),ylim=c(0.0,0.2),xaxt='n',
xlab='Time',ylab=expression(Connectance ~~ L/S^{2}))
"""
for t=2:length(seq_stitch)
    R"""
    points(jitter(rep($(t),$reps)),jitter($(conn_stitch[:,t])),pch='.',col=pal[2])
    """
end
R"""
boxplot($(conn_stitch),outline=FALSE,names=$(seq_stitch),xlab='',ylab='',
pars = list(boxwex = 0.2, staplewex = 0.5, outwex = 0.5,whisklty=1),add=TRUE)
points($(meanconn),ylim=c(0,0.03),pch=16)
lines($(meanconn),ylim=c(0,0.03),lwd=2)
lines(seq(0.001,3000),rep(mean($Pconn),3000),lty=3)
dev.off()
"""

mearlysp = mean(sprich[setdiff(collect(1:reps),earlygen),:],dims=1)
mearlygen = mean(sprich[earlygen,:],dims=1)

minrealizedG = Array{Float64}(undef,reps,length(seq2));
maxrealizedG = Array{Float64}(undef,reps,length(seq2));
for r=1:reps
    for t=1:length(seq2)
        minrealizedG[r,t] = minimum(realizedG[r,seq2[t],realizedG[r,seq2[t],:].>0]);
        maxrealizedG[r,t] = maximum(realizedG[r,seq2[t],realizedG[r,seq2[t],:].>0]);
    end
end

filename = "figures/yog/meanassembly_time.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1');
colpal = c(rep(pal[1],100),rep(pal[2],100))
pdf($namespace,height=10,width=8)
par(mfrow=c(3,1))
plot($seq,$(vec(mearlygen)),type='l',col=pal[1],ylim=c(0,150),log='x')
lines($seq,$(vec(mearlysp)),col=pal[2])
legend(2,150,legend=c('Early Gen', 'Early Spec'),col=c(pal[1],pal[2]),cex=0.8,pch=16)
plot($(seq[seq2]),$(vec(mean(rG[earlygen,:],dims=1))),ylim=c(0,2),type='l',col=pal[1],log='x')
points(jitter($(repeat(seq[seq2],outer=length(earlygen)))),$(vec(rG[earlygen,:])),pch='.',col=pal[1])
lines($(seq[seq2]),$(vec(mean(rG[earlysp,:],dims=1))),col=pal[2])
points(jitter($(repeat(seq[seq2],outer=length(earlysp)))),$(vec(rG[earlysp,:])),pch='.',col=pal[2])
plot(jitter($(repeat(seq[seq2],outer=length(earlygen)))),$(vec(minrealizedG[earlygen,:])),ylim=c(0,2),pch='.',col=pal[1],log='x')
points(jitter($(repeat(seq[seq2],outer=length(earlysp)))),$(vec(minrealizedG[earlysp,:])),pch='.',col=pal[2])
dev.off()
"""


filename = "figures/yog/specialistsvsgeneralists.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1');
colpal = c(rep(pal[1],100),rep(pal[2],100))
pdf($namespace,height=5,width=6)
par(mfrow=c(1,1))
plot($(seq[seq2]),$(mean(minrealizedG[earlygen,:],dims=1)),type='l',col='black',log='xy',ylim=c(0.5,4))
points($(seq[seq2]),$(mean(minrealizedG[earlygen,:],dims=1)),pch=21,col='black',bg=pal[1])
lines($(seq[seq2]),$(mean(minrealizedG[earlysp,:],dims=1)))
points($(seq[seq2]),$(mean(minrealizedG[earlysp,:],dims=1)),pch=21,col='black',bg=pal[2])

lines($(seq[seq2]),$(mean(maxrealizedG[earlygen,:],dims=1)))
points($(seq[seq2]),$(mean(maxrealizedG[earlygen,:],dims=1)),pch=21,col='black',bg=pal[1])
lines($(seq[seq2]),$(mean(maxrealizedG[earlysp,:],dims=1)))
points($(seq[seq2]),$(mean(maxrealizedG[earlysp,:],dims=1)),pch=21,col='black',bg=pal[2])

# boxplot($(minrealizedG[earlygen,:]),pch='.',col=pal[1],xlim=c(0,10),ylim=c(0,2))
# boxplot($(minrealizedG[earlysp,:]),pch='.',col=pal[2],add=TRUE)
# boxplot($(maxrealizedG[earlygen,:]),pch='.',col=pal[1],add=TRUE)
# boxplot($(maxrealizedG[earlysp,:]),pch='.',col=pal[2],add=TRUE)
dev.off()
"""


i=2;
#Visualize food web
r=earlygen[i];
tstep = seq[seq2[3]];
filename = "data/steadystate/int_m.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace int_m tp_m tind_m mp_m mind_m;
a_b,n_b,i_b,m_b,n_b0,sp_v,int_id = preamble_defs(int_m);
filename = "data/steadystate/cid.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace CID clock;

cid = findall(isodd,CID[:,tstep]);
deg,troph = structure(S,cid,sp_v,tind_m);
spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);
#Degree distribution
# degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
adjmatrix = tind_m[[1;spcid_ind],[1;spcid_ind]];
indmatrix = adjmatrix .- tp_m[[1;spcid_ind],[1;spcid_ind]];
dirmatrix = tp_m[[1;spcid_ind],[1;spcid_ind]];

filename = "figures/yog/foodweb_time2.pdf"
namespace = smartpath(filename);
R"""
library(igraph)
library(RColorBrewer)
pdf($namespace,width=8,height=5)
par(mfrow=c(1,2))
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(adjmatrix'));
basal_pos <- 1
trophic = as.numeric($([0;troph[1:size(adjmatrix)[1]-1]]));
#trophic = as.numeric($([0;paths[keepnodes[2:length(keepnodes)]]]));
keepnodes = c(1,which(trophic>0.9))"""; @rget keepnodes; keepnodes = Int64.(keepnodes);
R"""
#keepnodes = $keepnodes;
trophic2 = trophic[keepnodes];
coords <- cbind(runif(length(keepnodes)),trophic2);
coords[basal_pos,1] <- 0.5
fw_g = graph.adjacency($(adjmatrix[keepnodes,keepnodes]'))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)),main='Gen') 
#main=ecount(fw_g)/$(size(adjmatrix)[1])^2,
fw_ind <- graph.adjacency($(indmatrix[keepnodes,keepnodes]'));
#plot(fw_ind,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='red',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)),add=TRUE)
# dev.off()
"""
r=earlysp[i];
tstep = seq[seq2[3]];
filename = "data/steadystate/int_m.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace int_m tp_m tind_m mp_m mind_m;
a_b,n_b,i_b,m_b,n_b0,sp_v,int_id = preamble_defs(int_m);
filename = "data/steadystate/cid.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace CID clock;

cid = findall(isodd,CID[:,tstep]);
deg,troph = structure(S,cid,sp_v,tind_m);
spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);
#Degree distribution
# degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
adjmatrix = tind_m[[1;spcid_ind],[1;spcid_ind]];
indmatrix = adjmatrix .- tp_m[[1;spcid_ind],[1;spcid_ind]];
dirmatrix = tp_m[[1;spcid_ind],[1;spcid_ind]];
R"""
library(igraph)
library(RColorBrewer)
# pdf($namespace,width=6,height=5)
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(adjmatrix'));
basal_pos <- 1
trophic = as.numeric($([0;troph[1:size(adjmatrix)[1]-1]]));
#trophic = as.numeric($([0;paths[keepnodes[2:length(keepnodes)]]]));
keepnodes = c(1,which(trophic>0.9))"""; @rget keepnodes; keepnodes = Int64.(keepnodes);
R"""
#keepnodes = $keepnodes;
trophic2 = trophic[keepnodes];
coords <- cbind(runif(length(keepnodes)),trophic2);
coords[basal_pos,1] <- 0.5
fw_g = graph.adjacency($(adjmatrix[keepnodes,keepnodes]'))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)),main='Spec') 
#main=ecount(fw_g)/$(size(adjmatrix)[1])^2,
fw_ind <- graph.adjacency($(transpose(indmatrix[keepnodes,keepnodes])));
#plot(fw_ind,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='red',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)),add=TRUE)
dev.off()
"""

#WHAT IS GOING ON?
r=86
tstep = seq[seq2[3]];
filename = "data/steadystate/int_m.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace int_m tp_m tind_m mp_m mind_m;
a_b,n_b,i_b,m_b,n_b0,sp_v,int_id = preamble_defs(int_m);
filename = "data/steadystate/cid.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace CID clock;

cid = findall(isodd,CID[:,tstep]);
deg,troph = structure(S,cid,sp_v,tind_m);
spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);








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
# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/yog/trophic_degrees_time.pdf");
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
# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/yog/richconn.pdf");
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

# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Enigma/figures/yog/richconn2.pdf");
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
# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/yog/engcdf.pdf");
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
# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/yog/objcdf.pdf");
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

# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/yog/potcol.pdf");
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


