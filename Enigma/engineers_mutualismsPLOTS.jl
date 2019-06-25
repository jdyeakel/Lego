if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

reps = 50;
S = 200;
maxits = 4000;

nvec = collect(0.0:0.1:2.0);
lnvec = length(nvec);

lambdavec = collect(0:0.1:2.0);
llamb = length(lambdavec);

engineers = SharedArray{Int64}(lnvec,llamb,reps,maxits);
sprich = SharedArray{Int64}(lnvec,llamb,reps,maxits);
rich = SharedArray{Int64}(lnvec,llamb,reps,maxits);
clocks = SharedArray{Float64}(lnvec,llamb,reps,maxits);

mextrate = SharedArray{Float64}(lnvec,llamb,reps);
stdextrate = SharedArray{Float64}(lnvec,llamb,reps);

totalcolonizations = SharedArray{Float64}(lnvec,llamb,reps);
totalextinctions = SharedArray{Float64}(lnvec,llamb,reps);

mpersistance = SharedArray{Float64}(lnvec,llamb,reps);
stdpersistance = SharedArray{Float64}(lnvec,llamb,reps);



for v = 1:lnvec
    
    for w = 1:llamb
        
        filename = "data/engineers_mutualisms/sim_settings.jld";
        indices = [v,w];
        namespace = smartpath(filename,indices);
        @load namespace reps S maxits nvec athresh nthresh lambda SSprobs SOprobs OOprobs;
        
        @sync @distributed for r=1:reps
            
            filename = "data/engineers_mutualisms/int_m.jld";
            indices = [v,w,r];
            namespace = smartpath(filename,indices);
            @load namespace int_m tp_m tind_m mp_m mind_m;
            
            filename = "data/engineers_mutualisms/cid.jld";
            indices = [v,w,r];
            namespace = smartpath(filename,indices);
            @load namespace CID clock;
            
            a_b,
            n_b,
            i_b,
            m_b,
            n_b0,
            sp_v,
            int_id = preamble_defs(int_m);
            
            
            for t=1:maxits
                sprich[v,w,r,t] = sum(CID[1:S,t]);
                rich[v,w,r,t] = sum(CID[:,t]);
                #how many engineers?
                spcid = findall(isodd,CID[1:S,t]);
                engineers[v,w,r,t] = sum(sum(m_b[spcid,:],dims=2) .> 0);
            end

            clocks[v,w,r,:] = copy(clock);

            #Delta species and richness over time
            spdiff = diff(sprich[v,w,r,:]);
            rdiff = diff(rich[v,w,r,:]);
            dt = diff(clock);
            
            #Where Delta species > 0 (+1) is a colonization (position)
            colpos = findall(x->x>0,spdiff);
            #Where Delta species < 0 (-1 or less) is an extinction (position)
            extpos = findall(x->x<0,spdiff);
            #Extinction size (but make it positive)
            extinctions = spdiff[extpos].*-1;
            #Colonization size
            colonizations = spdiff[colpos];
            
            totalextinctions[v,w,r] = sum(extinctions);
            totalcolonizations[v,w,r] = sum(colonizations);
            
            tstepsincom = sum(CID[2:S,:],dims=2) ./ maxits;
            mpersistance[v,w,r] = mean(tstepsincom[vec(findall(!iszero,tstepsincom))]);
            stdpersistance[v,w,r] = std(tstepsincom[vec(findall(!iszero,tstepsincom))]);
            
            #Only loop if there are extinctions
            if length(extinctions) > 0

                #Number of extinctions / dt
                extrate = extinctions ./ dt[extpos];
                mextrate[v,w,r] = mean(extrate);
                stdextrate[v,w,r] = std(extrate);

                #Number of colonizations / dt
                colrate = colonizations ./ dt[colpos];
            else
                mextrate[v,w,r] = 0;
                stdextrate[v,w,r] = 0;
            end
            
        end
        println(string("v=",v,"/",lnvec,"; w=",w,"/",llamb))
    end
end
objects = rich .- sprich;



#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)
filename = "data/engineers_mutualisms/meanrates.jld";
namespace = smartpath(filename);
@save namespace reps lambdavec llamb sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations mpersistance stdpersistance


#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)
filename = "data/engineers_mutualisms/meanrates.jld";
namespace = smartpath(filename);
@load namespace reps lambdavec llamb sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations mpersistance stdpersistance

avalue = 0.01;
nvec_scaled = (avalue/10) .* nvec;
#mean extinction rate as a function of engineering and mutualisms
mextrate_surf = Array{Float64}(undef,lnvec,llamb);
# mcolrate_surf = Array{Float64}(undef,lnvec,llamb);
mpersist_surf = Array{Float64}(undef,lnvec,llamb);
for v = 1:lnvec
    for w = 1:llamb
        mextrate_surf[v,w] = mean(mextrate[v,w,:]);
        # mcolrate_surf[v,w] = mean(mcolrate[v,w,:]);
        mpersist_surf[v,w] = mean(mpersistance[v,w,:]);
    end
end


#x ~ row of z
#y ~ columns of z
filename = "figures/engmut/extrate_engmut.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pal = brewer.pal(11,'Spectral')
pdf($namespace,width=5,height=5)
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mextrate_surf)),col=pal,ylab='Frequency of mutualisms',xlab='Mean number of objects/species',nlevel=11)
# contour(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mextrate_surf)), add = TRUE, drawlabels = TRUE)

dev.off()
"""


#x ~ row of z
#y ~ columns of z
filename = "figures/engmut/persistence_engmut.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(9,'YlGnBu')
pdf($namespace,width=5,height=5)
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mpersist_surf)),col=pal,ylab='Frequency of mutualisms',xlab='Mean number of objects/species',nlevel=11)
# contour(x=$nvec_scaled,y=$lambdavec,z=$mpersist_surf, add = TRUE, drawlabels = TRUE)
dev.off()
"""


filename = "../manuscript/fig_engineers.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
library(fields)
pdf($namespace,width=5,height=5)
layout(matrix(c(1,2), 2, 1, byrow = TRUE), 
   widths=c(1,1), heights=c(0.5,0.5))
par(oma = c(0.5, 1, 1, 1), mar = c(3, 4, 0, 1))
pal = brewer.pal(11,'Spectral')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mextrate_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Frequency of service', line=2.0, cex.lab=1.)
title(xlab='Mean objects/species', line=1.5, cex.lab=1.)

pal = brewer.pal(9,'YlGnBu')
image.plot(y=$nvec_scaled,x=$lambdavec,z=$(transpose(mpersist_surf)),col=pal,ylab='',xlab='',nlevel=11,axes=FALSE)
axis(2,at=seq(0,0.002,by=0.001),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
axis(1,at=seq(0,2,by=0.5),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Frequency of service', line=2.0, cex.lab=1.)
title(xlab='Mean objects/species', line=1.5, cex.lab=1.)
dev.off()
"""
