#EXTINCTION RATES
if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end


filename = "data/engineers_unique/sim_settings.jld";
namespace = smartpath(filename);
@load namespace reps S maxits athresh nthresh lambdavec SSprobs SOprobs OOprobs;

llamb = length(lambdavec);
its = llamb*reps;

engineers = SharedArray{Int64}(its,maxits);
sprich = SharedArray{Int64}(its,maxits);
rich = SharedArray{Int64}(its,maxits);
clocks = SharedArray{Float64}(its,maxits);

mextrate = SharedArray{Float64}(its);
stdextrate = SharedArray{Float64}(its);

totalcolonizations = SharedArray{Float64}(its);
totalextinctions = SharedArray{Float64}(its);

mpersistance = SharedArray{Float64}(its);
stdpersistance = SharedArray{Float64}(its);

@sync @distributed for i = 0:(its - 1)
    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;

    ii = i + 1;
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    filename = "data/engineers_unique/int_m.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;


    filename = "data/engineers_unique/cid.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);


    #Calculate CDF of extinction cascade size

    #species richness over time
    #species + objects over time
    #engineers over time
    for t=1:maxits
        sprich[ii,t] = sum(CID[1:S,t]);
        rich[ii,t] = sum(CID[:,t]);
        #how many engineers?
        spcid = findall(isodd,CID[1:S,t]);
        engineers[ii,t] = sum(sum(m_b[spcid,:],dims=2) .> 0);
    end

    clocks[ii,:] = copy(clock);

    #Delta species and richness over time
    spdiff = diff(sprich[ii,:]);
    rdiff = diff(rich[ii,:]);

    dt = diff(clock);

    #Where Delta species > 0 (+1) is a colonization (position)
    colpos = findall(x->x>0,spdiff);
    #Where Delta species < 0 (-1 or less) is an extinction (position)
    extpos = findall(x->x<0,spdiff);
    #Extinction size (but make it positive)
    extinctions = spdiff[extpos].*-1;
    #Colonization size
    colonizations = spdiff[colpos];
    
    totalextinctions[ii] = sum(extinctions);
    totalcolonizations[ii] = sum(colonizations);
    
    tstepsincom = sum(CID[2:S,:],dims=2) ./ maxits;
    mpersistance[ii] = mean(tstepsincom[vec(findall(!iszero,tstepsincom))]);
    stdpersistance[ii] = std(tstepsincom[vec(findall(!iszero,tstepsincom))]);

    #Only loop if there are extinctions
    if length(extinctions) > 0

        #Number of extinctions / dt
        extrate = extinctions ./ dt[extpos];
        mextrate[ii] = mean(extrate);
        stdextrate[ii] = std(extrate);

        #Number of colonizations / dt
        colrate = colonizations ./ dt[colpos];

        # extratevec = collect(0:0.0001:maximum(extrate));
        # extratevec[ii,:] = collect(range(0,step=maximum(extrate)/lcdf,length=lcdf));

        # extcdf = Array{Int64}(undef,lcdf);
        # for j=1:lcdf
        #     extcdf[j] = length(findall(x->x<extratevec[ii,j],extrate));
        # end
        # 
        # EXTCDF[ii,:] = extcdf;
    else
        mextrate[ii] = 0;
        stdextrate[ii] = 0;
        # extratevec[ii,:] = repeat([0],inner=lcdf);
        # EXTCDF[ii,:] = repeat([0],inner=lcdf);
    end
    # if mod(r,1) == 0
    #     println("reps =",r)
    # end
end

#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)
filename = "data/engineers_unique/meanrates.jld";
namespace = smartpath(filename);
@save namespace reps lambdavec llamb sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations mpersistance stdpersistance;

#Load the CDF data
filename = "data/engineers_unique/meanrates.jld";
namespace = smartpath(filename);
@load namespace reps lambdavec llamb sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations;



its = reps*llamb;
objects = rich .- sprich;

#Mean steady state number of objects and engineers
mobj = vec(mean(objects[:,maxits-100:maxits],dims=2));
meng = vec(mean(engineers[:,maxits-100:maxits],dims=2));
mss = vec(mean(sprich[:,maxits-100:maxits],dims=2));
totaltime = clocks[:,4000];

#isolate diverse systems to minimize effects of s.s. richness
diverse = findall(x->x>120,mss);

#NOTE: this is just to test that reshape works the way it should!
mextarray_test = Array{Float64}(undef,reps,llamb)
for i = 0:(its - 1)
    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;
    
    ii=i+1;
    mextarray_test[b,a] = mextrate[ii];
end
        

#Mean extinction rate vs. SS engineers
filename = "figures/eng/mean_extrate_unique.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,width=8,height=4)
par(mfrow=c(1,2))
plot($(meng[diverse]),$(mextrate[diverse]),pch='.',ylab='Mean extinction rate',xlab='Mean num. engineers')
plot($(mss[diverse]),$(mextrate[diverse]),pch='.',ylab='Mean extinction rate',xlab='S.S. species richness')
dev.off()
"""

#Mean extinction rate vs. SS engineers
filename = "figures/eng/mean_box_extrate_unique.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,width=8,height=4)
par(mfrow=c(1,2))
plot($(meng[diverse]),$(mextrate[diverse]),pch='.',ylab='Mean extinction rate',xlab='Mean num. engineers')
boxplot($mextarray,names = $lambdavec,xlab='Mean num. objects/species',ylab='Mean extinction rate',outline=F,col='gray')
dev.off()
"""

###############################
#NOTE EXTINCTION RATE FIGURES
##############################

mextarray = reshape(Array(mextrate),reps,llamb);
filename = "figures/eng/boxplot_extrate_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($mextarray,names = $lambdavec,xlab='Mean number of objects/species',ylab='Extinciton rate mean',outline=F,col='gray')
dev.off()
"""


mssarray = reshape(Array(mss),reps,llamb);
filename = "figures/eng/boxplot_mss_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($mssarray,names = $lambdavec,xlab='Mean number of objects/species',ylab='Species richness',outline=F,col='gray')
dev.off()
"""

mextarray = reshape(Array(mextrate),reps,llamb);
mssarray = reshape(Array(mss),reps,llamb);
filename = "figures/eng/boxplot_extratepercapita_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($mextarray / $mssarray,names = $lambdavec,xlab='Mean number of objects/species',ylab='Extinction rate (per capita)',outline=F,col='gray')
dev.off()
"""


totextinctions = reshape(Array(totalextinctions),reps,llamb);
tottime = reshape(Array(totaltime),reps,llamb);
mssarray = reshape(Array(mss),reps,llamb);
filename = "figures/eng/boxplot_totext_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($totextinctions / $tottime,names = $lambdavec,xlab='Mean number of objects/species',ylab='Total extinctions / Total time',outline=F,col='gray')
dev.off()
"""

filename = "figures/eng/boxplot_totextmss_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot(($totextinctions / ($tottime * $mssarray)),names = $lambdavec,xlab='Mean number of objects/species',ylab='Total extinctions / (Total time * SS)',outline=F,col='gray')
dev.off()
"""



mextarray = reshape(Array(mextrate),reps,llamb);
stdextarray = reshape(Array(stdextrate),reps,llamb);
filename = "figures/eng/boxplot_std_extrate_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($stdextarray,names = $lambdavec,xlab='Mean number of objects/species',ylab='Extinciton rate SD',outline=F,col='gray')
dev.off()
"""

mextarray = reshape(Array(mextrate),reps,llamb);
stdextarray = reshape(Array(stdextrate),reps,llamb);
filename = "figures/eng/boxplot_cv_extrate_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($stdextarray / $mextarray,names = $lambdavec,xlab='Mean number of objects/species',ylab='Extinciton rate CV',outline=F,col='gray')
dev.off()
"""

#Bulk number of extinctions
extarray = reshape(Array(totalextinctions),reps,llamb);
colarray = reshape(Array(totalcolonizations),reps,llamb);
filename = "figures/eng/boxplot_col_ext_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=10)
par(mfrow=c(2,1))
boxplot($colarray/$(reshape(Array(mss),reps,llamb)),names = $lambdavec,xlab='Mean number of objects/species',ylab='Total colonizations',outline=F,col='gray')
boxplot($extarray/$(reshape(Array(mss),reps,llamb)),names = $lambdavec,xlab='Mean number of objects/species',ylab='Total Extinctions',outline=F,col='gray')
dev.off()
"""

filename = "figures/eng/boxplot_persistance_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($(reshape(Array(mpersistance),reps,llamb)),names = $lambdavec,xlab='Mean number of objects/species',ylab='Mean persistance (%)',outline=F,col='gray',ylim=c(0.5,1))
dev.off()
"""


filename = "figures/eng/boxplot_sdpersistance_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($(reshape(Array(stdpersistance),reps,llamb)),names = $lambdavec,xlab='Mean number of objects/species',ylab='Persistance SD (%)',outline=F,col='gray',ylim=c(0,0.5))
dev.off()
"""

filename = "figures/eng/boxplot_cvpersistance_unique.pdf"
namespace = smartpath(filename)
R"""
pdf($namespace,width=6,height=5)
boxplot($(reshape(Array(stdpersistance),reps,llamb))/$(reshape(Array(mpersistance),reps,llamb)),names = $lambdavec,xlab='Mean number of objects/species',ylab='Persistance CV',outline=F,col='gray',ylim=c(0,1))
dev.off()
"""





















#POTENTIAL COLONIZERS


filename = "data/engineers_unique/sim_settings.jld";
namespace = smartpath(filename);
@load namespace reps S maxits athresh nthresh lambdavec SSprobs SOprobs OOprobs;

llamb = length(lambdavec);
its = llamb*reps;


#POTENTIAL COLONIZERS
sprich = SharedArray{Int64}(llamb,reps,maxits);
rich = SharedArray{Int64}(llamb,reps,maxits);
pc = SharedArray{Int64}(llamb,reps,maxits);
@sync @distributed for i = 0:(its - 1)

    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;

    ii = i + 1;
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    filename = "data/engineers_unique/int_m.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;


    filename = "data/engineers_unique/cid.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;

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
        sprich[a,b,t] = sum(CID[1:S,t]);
        rich[a,b,t] = sum(CID[:,t]);
        pcid, lpc = potcol(sp_v,int_id,cid,a_b,n_b0,athresh,nthresh);
        pc[a,b,t] = lpc;
    end

end


filename = "data/engineers_unique/potcol.jld";
namespace = smartpath(filename);
@save namespace sprich rich pc;
#
# save(namespace,
# "sprich", sprich,
# "rich", rich,
# "pc", pc
# );
#
# dpc = load(namespace);
# sprich = dpc["sprich"],
# rich = dpc["rich"],
# pc = dpc["pc"],
# );

#NOTE: This needs changed to work =)

loweng = findall(x->x==0.5,lambdavec)[1];
medeng = findall(x->x==1.0,lambdavec)[1];
higheng = findall(x->x==2.0,lambdavec)[1];

mpc = vec(mean(pc[loweng,:,:],dims=1));
sdpc = vec(std(pc[loweng,:,:],dims=1));
propss = vec(mean(sprich[loweng,:,:],dims=1)) ./ mean(sprich[loweng,:,maxits-100:maxits]);
ns = mean(sprich[loweng,:,maxits-100:maxits]);
filename = "figures/eng/potcol22_unique.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pdf($namespace,width=6,height=5)
pal = brewer.pal(3,'Set1');
plot($propss,$mpc/$ns,type='l',col=pal[1],lwd=2,xlab='Proportion filled',ylab='Available niche space',ylim=c(0,0.5))
polygon(x=c($propss,rev($propss)),y=c(($mpc-$sdpc)/$ns,rev(($mpc+$sdpc)/$ns)),col=paste(pal[1],50,sep=''),border=NA)
lines($propss,$mpc/$ns,col=pal[1],lwd=2)
"""
mpc = vec(mean(pc[medeng,:,:],dims=1));
sdpc = vec(std(pc[medeng,:,:],dims=1));
propss = vec(mean(sprich[medeng,:,:],dims=1)) ./ mean(sprich[medeng,:,maxits-100:maxits]);
ns = mean(sprich[medeng,:,maxits-100:maxits]);
R"""
polygon(x=c($propss,rev($propss)),y=c(($mpc-$sdpc)/$ns,rev(($mpc+$sdpc)/$ns)),col=paste(pal[2],50,sep=''),border=NA)
lines($propss,$mpc/$ns,col=pal[2],lwd=2)
# dev.off()
"""
mpc = vec(mean(pc[higheng,:,:],dims=1));
sdpc = vec(std(pc[higheng,:,:],dims=1));
propss = vec(mean(sprich[higheng,:,:],dims=1)) ./ mean(sprich[higheng,:,maxits-100:maxits]);
ns = mean(sprich[higheng,:,maxits-100:maxits]);
R"""
polygon(x=c($propss,rev($propss)),y=c(($mpc-$sdpc)/$ns,rev(($mpc+$sdpc)/$ns)),col=paste(pal[3],50,sep=''),border=NA)
lines($propss,$mpc/$ns,col=pal[3],lwd=2)
dev.off()
"""

#BOTH MEASURES NORMALIZED TO STEADY STATE
mpc = vec(mean(pc[loweng,:,:],dims=1));
sdpc = vec(std(pc[loweng,:,:],dims=1));
propss = vec(mean(sprich[loweng,:,:],dims=1)) ./ mean(sprich[loweng,:,maxits-100:maxits]);
ns = mean(sprich[loweng,:,maxits-100:maxits]);
filename = "figures/eng/potcol32_unique.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pdf($namespace,width=6,height=5)
pal = brewer.pal(3,'Set1');
plot($propss,$mpc/$ns,type='l',col=pal[1],lwd=2,xlab='Proportion filled',ylab='Available niche space',ylim=c(0,0.5))
"""
mpc = vec(mean(pc[medeng,:,:],dims=1));
sdpc = vec(std(pc[medeng,:,:],dims=1));
propss = vec(mean(sprich[medeng,:,:],dims=1)) ./ mean(sprich[medeng,:,maxits-100:maxits]);
ns = mean(sprich[medeng,:,maxits-100:maxits]);
R"""
lines($propss,$mpc/$ns,col=pal[2],lwd=2)
# dev.off()
"""
mpc = vec(mean(pc[higheng,:,:],dims=1));
sdpc = vec(std(pc[higheng,:,:],dims=1));
propss = vec(mean(sprich[higheng,:,:],dims=1)) ./ mean(sprich[higheng,:,maxits-100:maxits]);
ns = mean(sprich[higheng,:,maxits-100:maxits]);
R"""
lines($propss,$mpc/$ns,col=pal[3],lwd=2)
dev.off()
"""


#The 'diverse' code currently won't work as written
mss = vec(mean(sprich[:,:,maxits-100:maxits],dims=2));
diverse = findall(x->x>120,mss);


#UN-NORMALIZED
mpc = vec(mean(pc[loweng,:,:],dims=1));
sdpc = vec(std(pc[loweng,:,:],dims=1));
propss = vec(mean(sprich[loweng,:,:],dims=1)) ./ mean(sprich[loweng,:,maxits-100:maxits]);
ns = mean(sprich[loweng,:,maxits-100:maxits]);
filename = "figures/eng/potcol32_unnorm_unique.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pdf($namespace,width=6,height=5)
pal = brewer.pal(3,'Set1');
plot($propss,$mpc,type='l',col=pal[1],lwd=2,xlab='Proportion filled',ylab='Available niche space',ylim=c(0,0.5))
"""
mpc = vec(mean(pc[medeng,:,:],dims=1));
sdpc = vec(std(pc[medeng,:,:],dims=1));
propss = vec(mean(sprich[medeng,:,:],dims=1)) ./ mean(sprich[medeng,:,maxits-100:maxits]);
ns = mean(sprich[medeng,:,maxits-100:maxits]);
R"""
lines($propss,$mpc,col=pal[2],lwd=2)
# dev.off()
"""
mpc = vec(mean(pc[higheng,:,:],dims=1));
sdpc = vec(std(pc[higheng,:,:],dims=1));
propss = vec(mean(sprich[higheng,:,:],dims=1)) ./ mean(sprich[higheng,:,maxits-100:maxits]);
ns = mean(sprich[higheng,:,maxits-100:maxits]);
R"""
lines($propss,$mpc,col=pal[3],lwd=2)
dev.off()
"""
