#EXTINCTION RATES
if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end


filename = "data/engineers2/sim_settings.jld";
namespace = smartpath(filename);
@load namespace reps S maxits athresh nthresh lambdavec SSprobs SOprobs OOprobs;

llamb = length(lambdavec);
its = llamb*reps;

numcol = SharedArray{Int64}(llamb,reps);
numprim = SharedArray{Int64}(llamb,reps);
numsec = SharedArray{Int64}(llamb,reps);
numobext = SharedArray{Int64}(llamb,reps);
clocks = SharedArray{Float64}(llamb,reps);
mss = SharedArray{Float64}(llamb,reps);
mprimextrate = SharedArray{Float64}(llamb,reps);
msecextrate = SharedArray{Float64}(llamb,reps);
mcolrate = SharedArray{Float64}(llamb,reps);
mobextrate = SharedArray{Float64}(llamb,reps);

@sync @distributed for i = 0:(its - 1)
    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;

    ii = i + 1;
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    filename = "data/engineers2/int_m.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;


    filename = "data/engineers2/cid.jld";
    indices = [a,b];
    namespace = smartpath(filename,indices);
    @load namespace CID clock events;

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    colpos = findall(x->x==1,events);
    primextpos = findall(x->x==1,events);
    secextpos = findall(x->x==2,events);
    obextpos = findall(x->x==3,events);
    #number of primary extinctions
    numcol[a,b] = length(colpos);
    numprim[a,b] = length(primextpos);
    numsec[a,b] = length(secextpos);
    numobext[a,b] = length(obextpos);
    clocks[a,b] = clock[maxits];
    mss[a,b] = mean(sum(CID[:,maxits-100:maxits]));
    
    mprimextrate[a,b] = mean(1 ./ clock[primextpos]);
    msecextrate[a,b] = mean(1 ./ clock[secextpos]);
    mcolrate[a,b] = mean(1 ./ clock[colpos]);
    mobextrate[a,b] = mean(1 ./ clock[obextpos]);
    
end

filename = "figures/eng2/fig_extinctionrates.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,width=5,height=4)
plot($lambdavec,$(mean(mprimextrate,dims=2)),pch=16,col=pal[1],ylim=c(0.02,0.06))
points($lambdavec,$(mean(msecextrate,dims=2)),pch=16,col=pal[2])
dev.off()
"""



filename = "figures/eng2/fig_extinctionratio.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,width=5,height=4)
plot($lambdavec,$(mean(mprimextrate ./ msecextrate,dims=2)),pch=16,col=pal[1])
dev.off()
"""


filename = "figures/eng2/fig_numext.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,width=5,height=4)
plot($lambdavec,$(mean(numprim ./ clocks,dims=2)),pch=16,col=pal[1],ylim=c(0,25),xlab="Number objects/species",ylab="Total extinctions/Total time")
points($lambdavec,$(mean(numsec ./ clocks,dims=2)),pch=16,col=pal[2])
dev.off()
"""



filename = "figures/eng2/fig_numextratio.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,width=5,height=4)
plot($lambdavec,$(mean((numprim ./ clocks) ./ (numsec ./ clocks),dims=2)),pch=16,col=pal[1],xlab="Number objects/species",ylab="1°/2° extinctions")
dev.off()
"""



filename = "figures/eng2/fig_numextratio_box.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal=brewer.pal(3,'Set1')
pdf($namespace,width=5,height=4)
boxplot($(transpose((numprim ./ clocks) ./ (numsec ./ clocks))),names=$lambdavec,pch=16,col=pal[1],xlab="Number objects/species",ylab="1°/2° extinctions",ylim=c(1,5))
dev.off()
"""

