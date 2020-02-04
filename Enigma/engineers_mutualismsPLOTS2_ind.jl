if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

filename = "data/engineers_mutualisms2/sim_settings.jld";
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
# numprim = SharedArray{Int64}(lnvec,llamb,reps);
# numsec = SharedArray{Int64}(lnvec,llamb,reps);
# numobext = SharedArray{Int64}(lnvec,llamb,reps);
# clocks = SharedArray{Float64}(lnvec,llamb,reps);
# mss = SharedArray{Float64}(lnvec,llamb,reps);
# mprimextrate = SharedArray{Float64}(lnvec,llamb,reps);
# msecextrate = SharedArray{Float64}(lnvec,llamb,reps);
# mcolrate = SharedArray{Float64}(lnvec,llamb,reps);
# mobextrate = SharedArray{Float64}(lnvec,llamb,reps);
# mpersistance = SharedArray{Float64}(lnvec,llamb,reps);
# propcol = SharedArray{Float64}(lnvec,llamb,reps);
# 
# numcol_unique = SharedArray{Int64}(lnvec,llamb,reps);
# numprim_unique = SharedArray{Int64}(lnvec,llamb,reps);
# numsec_unique = SharedArray{Int64}(lnvec,llamb,reps);
# numobext_unique = SharedArray{Int64}(lnvec,llamb,reps);
# clocks_unique = SharedArray{Float64}(lnvec,llamb,reps);
# mss_unique = SharedArray{Float64}(lnvec,llamb,reps);
# mprimextrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
# msecextrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
# mcolrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
# mobextrate_unique = SharedArray{Float64}(lnvec,llamb,reps);
# mpersistance_unique = SharedArray{Float64}(lnvec,llamb,reps);
# propcol_unique = SharedArray{Float64}(lnvec,llamb,reps);
# 
# lvec = copy(lnvec);
# paramvec = Array{Int64}(undef,lvec*lvec*reps,3);
# paramvec[:,1] = repeat(collect(1:lvec),inner=lvec*reps);
# paramvec[:,2] = repeat(collect(1:lvec),inner=reps,outer=lvec);
# paramvec[:,3] = repeat(collect(1:reps),outer=lvec*lvec);


#Look at:
# v = 1;
# w = 1, 11, 41

S = 200;

engvec = [1, 11, 41];
lengvec = length(engvec);

speciespersist = SharedArray{Float64}(lengvec,reps,S-1);
pindegree = SharedArray{Float64}(lengvec,reps,S-1);
poutdegree = SharedArray{Float64}(lengvec,reps,S-1);
tpindegree = SharedArray{Float64}(lengvec,reps,S-1);
tpoutdegree = SharedArray{Float64}(lengvec,reps,S-1);
mpindegree = SharedArray{Float64}(lengvec,reps,S-1);
mpoutdegree = SharedArray{Float64}(lengvec,reps,S-1);

routdegree = SharedArray{Float64}(lengvec,reps,S-1);
rindegree = SharedArray{Float64}(lengvec,reps,S-1);
troutdegree = SharedArray{Float64}(lengvec,reps,S-1);
trindegree = SharedArray{Float64}(lengvec,reps,S-1);
mroutdegree = SharedArray{Float64}(lengvec,reps,S-1);
mrindegree = SharedArray{Float64}(lengvec,reps,S-1);


for i=1:length(engvec)
    w = engvec[i];
    v = 1;
    for r=1:reps
        
        filename = "data/engineers_mutualisms2/int_m.jld";
        indices = [v,w,r];
        namespace = smartpath(filename,indices);
        @load namespace int_m tp_m tind_m mp_m mind_m;
        
        filename = "data/engineers_mutualisms2/cid.jld";
        indices = [v,w,r];
        namespace = smartpath(filename,indices);
        @load namespace CID clock events;
        
        colpos = findall(x->x==1,events);
        primextpos = findall(x->x==1,events);
        secextpos = findall(x->x==2,events);
        obextpos = findall(x->x==3,events);
        
        #persistence of each species
        tstepsincom = sum(CID[2:S,:],dims=2) ./ maxits;
        speciespersist[i,r,1:(S-1)] = tstepsincom;
        
        #POTENTIAL out-degree of each species
        poutdegree[i,r,1:(S-1)] = sum(tp_m[2:S,2:S] .+ mp_m[2:S,2:S],dims=1);
        #POTENTIAL in-degree of each species
        pindegree[i,r,1:(S-1)] = sum(tp_m[2:S,2:S] .+ mp_m[2:S,2:S],dims=2);
        
        tpoutdegree[i,r,1:(S-1)] = sum(tp_m[2:S,2:S],dims=1);
        #POTENTIAL in-degree of each species
        tpindegree[i,r,1:(S-1)] = sum(tp_m[2:S,2:S],dims=2);
        
        mpoutdegree[i,r,1:(S-1)] = sum(mp_m[2:S,2:S],dims=1);
        #POTENTIAL in-degree of each species
        mpindegree[i,r,1:(S-1)] = sum(mp_m[2:S,2:S],dims=2);
        
        
        
        #temporal analysis
        rod = zeros(Float64,S,maxits);
        rid = zeros(Float64,S,maxits);
        trod = zeros(Float64,S,maxits);
        trid = zeros(Float64,S,maxits);
        mrod = zeros(Float64,S,maxits);
        mrid = zeros(Float64,S,maxits);
        for t=1:maxits
            cid = findall(isodd,CID[:,t]);
            spcid = intersect(collect(2:S),cid);
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
        
        routdegree[i,r,:] = vec(meanfinite(rod,2))[2:S];
        rindegree[i,r,:] = vec(meanfinite(rid,2))[2:S];
        
        troutdegree[i,r,:] = vec(meanfinite(trod,2))[2:S];
        trindegree[i,r,:] = vec(meanfinite(trid,2))[2:S];
        
        mroutdegree[i,r,:] = vec(meanfinite(mrod,2))[2:S];
        mrindegree[i,r,:] = vec(meanfinite(mrid,2))[2:S];
        
    end
end



namespace = "$(homedir())/2014_Lego/Enigma/figures/yog/fig_indeng1.pdf";
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=10,height=6)
layout(matrix(c(1,3,2,4),2,2,byrow=TRUE))
par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
"""

eng = 1;
m1pos_in = findall(x-> x==1,vec(trindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(trindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(trindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(trindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(trindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(trindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(trindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(trindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(troutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(troutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(troutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(troutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(troutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(troutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(troutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(troutdegree[eng,:,:]));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
# dev.off()
"""

eng = 1;
m1pos_in = findall(x-> x==1,vec(tpindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(tpindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(tpindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(tpindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(tpindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(tpindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(tpindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(tpindegree[eng,:,:]));
m9pos_in = findall(x->x==9,vec(tpindegree[eng,:,:]));
m10pos_in = findall(x->x==10,vec(tpindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(tpoutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(tpoutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(tpoutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(tpoutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(tpoutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(tpoutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(tpoutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(tpoutdegree[eng,:,:]));
m9pos_out = findall(x->x==9,vec(tpoutdegree[eng,:,:]));
m10pos_out = findall(x->x==10,vec(tpoutdegree[eng,:,:]));
# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]),$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in]),$(speciespersist[eng,:,:][m9pos_in]),$(speciespersist[eng,:,:][m10pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))),rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in]))),rep(9,length($(speciespersist[eng,:,:][m9pos_in]))),rep(10,length($(speciespersist[eng,:,:][m10pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]),$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out]),$(speciespersist[eng,:,:][m9pos_out]),$(speciespersist[eng,:,:][m10pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))),rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out]))),rep(9,length($(speciespersist[eng,:,:][m9pos_out]))),rep(10,length($(speciespersist[eng,:,:][m10pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
dev.off()
"""


##############
# ENG2
#############

namespace = "$(homedir())/2014_Lego/Enigma/figures/yog/fig_indeng2.pdf";
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=10,height=6)
layout(matrix(c(1,3,2,4),2,2,byrow=TRUE))
par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
"""

eng = 2;
m1pos_in = findall(x-> x==1,vec(trindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(trindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(trindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(trindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(trindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(trindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(trindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(trindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(troutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(troutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(troutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(troutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(troutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(troutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(troutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(troutdegree[eng,:,:]));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
# dev.off()
"""

eng = 2;
m1pos_in = findall(x-> x==1,vec(tpindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(tpindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(tpindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(tpindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(tpindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(tpindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(tpindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(tpindegree[eng,:,:]));
m9pos_in = findall(x->x==9,vec(tpindegree[eng,:,:]));
m10pos_in = findall(x->x==10,vec(tpindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(tpoutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(tpoutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(tpoutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(tpoutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(tpoutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(tpoutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(tpoutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(tpoutdegree[eng,:,:]));
m9pos_out = findall(x->x==9,vec(tpoutdegree[eng,:,:]));
m10pos_out = findall(x->x==10,vec(tpoutdegree[eng,:,:]));
# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]),$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in]),$(speciespersist[eng,:,:][m9pos_in]),$(speciespersist[eng,:,:][m10pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))),rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in]))),rep(9,length($(speciespersist[eng,:,:][m9pos_in]))),rep(10,length($(speciespersist[eng,:,:][m10pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]),$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out]),$(speciespersist[eng,:,:][m9pos_out]),$(speciespersist[eng,:,:][m10pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))),rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out]))),rep(9,length($(speciespersist[eng,:,:][m9pos_out]))),rep(10,length($(speciespersist[eng,:,:][m10pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
dev.off()
"""


##############
# ENG3
#############

namespace = "$(homedir())/2014_Lego/Enigma/figures/yog/fig_indeng3.pdf";
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=10,height=6)
layout(matrix(c(1,3,2,4),2,2,byrow=TRUE))
par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
"""

eng = 3;
m1pos_in = findall(x-> x==1,vec(trindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(trindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(trindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(trindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(trindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(trindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(trindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(trindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(troutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(troutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(troutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(troutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(troutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(troutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(troutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(troutdegree[eng,:,:]));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
# dev.off()
"""

eng = 3;
m1pos_in = findall(x-> x==1,vec(tpindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(tpindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(tpindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(tpindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(tpindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(tpindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(tpindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(tpindegree[eng,:,:]));
m9pos_in = findall(x->x==9,vec(tpindegree[eng,:,:]));
m10pos_in = findall(x->x==10,vec(tpindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(tpoutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(tpoutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(tpoutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(tpoutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(tpoutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(tpoutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(tpoutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(tpoutdegree[eng,:,:]));
m9pos_out = findall(x->x==9,vec(tpoutdegree[eng,:,:]));
m10pos_out = findall(x->x==10,vec(tpoutdegree[eng,:,:]));
# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]),$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in]),$(speciespersist[eng,:,:][m9pos_in]),$(speciespersist[eng,:,:][m10pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))),rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in]))),rep(9,length($(speciespersist[eng,:,:][m9pos_in]))),rep(10,length($(speciespersist[eng,:,:][m10pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]),$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out]),$(speciespersist[eng,:,:][m9pos_out]),$(speciespersist[eng,:,:][m10pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))),rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out]))),rep(9,length($(speciespersist[eng,:,:][m9pos_out]))),rep(10,length($(speciespersist[eng,:,:][m10pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Pot. Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
dev.off()
"""




#######################################
#COMBINED PLOT (JUST REALIZED)
#######################################


namespace = "$(homedir())/2014_Lego/Enigma/figures/yog/fig_indeng_combined.pdf";
R"""
library(RColorBrewer)
pal = brewer.pal(3,'Set1')
pdf($namespace,width=6,height=4)
layout(matrix(c(1,3,5,2,4,6),2,3,byrow=TRUE))
par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
"""

eng = 1;
m1pos_in = findall(x-> x==1,vec(trindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(trindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(trindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(trindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(trindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(trindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(trindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(trindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(troutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(troutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(troutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(troutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(troutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(troutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(troutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(troutdegree[eng,:,:]));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
# dev.off()
"""


eng = 2;
m1pos_in = findall(x-> x==1,vec(trindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(trindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(trindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(trindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(trindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(trindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(trindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(trindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(troutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(troutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(troutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(troutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(troutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(troutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(troutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(troutdegree[eng,:,:]));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
# dev.off()
"""


eng = 3;
m1pos_in = findall(x-> x==1,vec(trindegree[eng,:,:]));
m2pos_in = findall(x->x==2,vec(trindegree[eng,:,:]));
m3pos_in = findall(x->x==3,vec(trindegree[eng,:,:]));
m4pos_in = findall(x->x==4,vec(trindegree[eng,:,:]));
m5pos_in = findall(x->x==5,vec(trindegree[eng,:,:]));
m6pos_in = findall(x->x==6,vec(trindegree[eng,:,:]));
m7pos_in = findall(x->x==7,vec(trindegree[eng,:,:]));
m8pos_in = findall(x->x==8,vec(trindegree[eng,:,:]));
m1pos_out = findall(x-> x==1,vec(troutdegree[eng,:,:]));
m2pos_out = findall(x->x==2,vec(troutdegree[eng,:,:]));
m3pos_out = findall(x->x==3,vec(troutdegree[eng,:,:]));
m4pos_out = findall(x->x==4,vec(troutdegree[eng,:,:]));
m5pos_out = findall(x->x==5,vec(troutdegree[eng,:,:]));
m6pos_out = findall(x->x==6,vec(troutdegree[eng,:,:]));
m7pos_out = findall(x->x==7,vec(troutdegree[eng,:,:]));
m8pos_out = findall(x->x==8,vec(troutdegree[eng,:,:]));

# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/fig_persistdegree_boxall.pdf";
R"""
# pdf($namespace,width=6,height=8)
# layout(matrix(c(2,2),2,1,byrow=TRUE))
# par(list(oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 1)))
data = c($(speciespersist[eng,:,:][m1pos_in]),$(speciespersist[eng,:,:][m2pos_in]),$(speciespersist[eng,:,:][m3pos_in]),$(speciespersist[eng,:,:][m4pos_in]))
# ,$(speciespersist[eng,:,:][m5pos_in]),$(speciespersist[eng,:,:][m6pos_in]),$(speciespersist[eng,:,:][m7pos_in]),$(speciespersist[eng,:,:][m8pos_in])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_in]))),rep(2,length($(speciespersist[eng,:,:][m2pos_in]))),rep(3,length($(speciespersist[eng,:,:][m3pos_in]))),rep(4,length($(speciespersist[eng,:,:][m4pos_in]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_in]))),rep(6,length($(speciespersist[eng,:,:][m6pos_in]))),rep(7,length($(speciespersist[eng,:,:][m7pos_in]))),rep(8,length($(speciespersist[eng,:,:][m8pos_in])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic in-degree (prey)',ylab='Persistence',ylim=c(0,1),las=1)

data = c($(speciespersist[eng,:,:][m1pos_out]),$(speciespersist[eng,:,:][m2pos_out]),$(speciespersist[eng,:,:][m3pos_out]),$(speciespersist[eng,:,:][m4pos_out]))
# ,$(speciespersist[eng,:,:][m5pos_out]),$(speciespersist[eng,:,:][m6pos_out]),$(speciespersist[eng,:,:][m7pos_out]),$(speciespersist[eng,:,:][m8pos_out])
f <-c(rep(1,length($(speciespersist[eng,:,:][m1pos_out]))),rep(2,length($(speciespersist[eng,:,:][m2pos_out]))),rep(3,length($(speciespersist[eng,:,:][m3pos_out]))),rep(4,length($(speciespersist[eng,:,:][m4pos_out]))))
# ,rep(5,length($(speciespersist[eng,:,:][m5pos_out]))),rep(6,length($(speciespersist[eng,:,:][m6pos_out]))),rep(7,length($(speciespersist[eng,:,:][m7pos_out]))),rep(8,length($(speciespersist[eng,:,:][m8pos_out])))
df = data.frame(x=data,f1 = f)
df$f1<-factor(df$f1)
boxplot(df$x ~ df$f1,col=pal[1],boxwex=0.25,xlab='Trophic out-degree (predators)',ylab='Persistence',ylim=c(0,1),las=1)
dev.off()
"""
