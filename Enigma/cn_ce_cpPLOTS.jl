if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

reps = 100;
S = 200;
maxits = 4000;

cnvec = collect(0:0.1:3.0);
cevec = copy(cnvec);
cp = 1.0;
lvec = length(cnvec);

paramvec = Array{Int64}(undef,lvec*lvec*reps,3);
paramvec[:,1] = repeat(collect(1:lvec),inner=lvec*reps);
paramvec[:,2] = repeat(collect(1:lvec),inner=reps,outer=lvec);
paramvec[:,3] = repeat(collect(1:reps),outer=lvec*lvec);


engineers = SharedArray{Int64}(lvec,lvec,reps,maxits);
sprich = SharedArray{Int64}(lvec,lvec,reps,maxits);
rich = SharedArray{Int64}(lvec,lvec,reps,maxits);
clocks = SharedArray{Float64}(lvec,lvec,reps,maxits);

mextrate = SharedArray{Float64}(lvec,lvec,reps);
stdextrate = SharedArray{Float64}(lvec,lvec,reps);

totalcolonizations = SharedArray{Float64}(lvec,lvec,reps);
totalextinctions = SharedArray{Float64}(lvec,lvec,reps);

mpersistance = SharedArray{Float64}(lvec,lvec,reps);
stdpersistance = SharedArray{Float64}(lvec,lvec,reps);

@sync @distributed for ii=1:(lvec*lvec*reps)
    
    v = paramvec[ii,1];
    w = paramvec[ii,2];
    r = paramvec[ii,3];
        
    filename = "data/cn_ce_cp/sim_settings.jld";
    indices = [v,w];
    namespace = smartpath(filename,indices);
    @load namespace reps S maxits cnvec athresh nthresh lambda SSprobs SOprobs OOprobs;
        
    
            
    filename = "data/cn_ce_cp/int_m.jld";
    indices = [v,w,r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    
    filename = "data/cn_ce_cp/cid.jld";
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
objects = rich .- sprich;



#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)
filename = "data/cn_ce_cp/meanrates.jld";
namespace = smartpath(filename);
@save namespace reps lvec sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations mpersistance stdpersistance


#SO WE DON'T HAVE TO RUN THE ABOVE ANALYSIS EVERY TIME (takes long time)
filename = "data/cn_ce_cp/meanrates.jld";
namespace = smartpath(filename);
@load namespace reps lvec sprich rich clocks engineers maxits mextrate stdextrate totalextinctions totalcolonizations mpersistance stdpersistance


#mean extinction rate as a function of engineering and mutualisms
mextrate_surf = Array{Float64}(undef,lvec,lvec);
# mcolrate_surf = Array{Float64}(undef,lvec,lvec);
mpersist_surf = Array{Float64}(undef,lvec,lvec);
for v = 1:lvec
    for w = 1:lvec
        mextrate_surf[v,w] = mean(mextrate[v,w,:]);
        # mcolrate_surf[v,w] = mean(mcolrate[v,w,:]);
        mpersist_surf[v,w] = mean(mpersistance[v,w,:]);
    end
end


#x ~ row of z
#y ~ columns of z
filename = "figures/extrate_cncecp.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(11,'Spectral')
pdf($namespace,width=5,height=5)
image(x=$cnvec,y=$cevec,z=$mextrate_surf,col=pal,xlab='cn',ylab='ce')
contour(x=$cnvec,y=$cevec,z=$mextrate_surf, add = TRUE, drawlabels = TRUE)
dev.off()
"""


#x ~ row of z
#y ~ columns of z
filename = "figures/persistence_cncecp.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(9,'YlGnBu')
pdf($namespace,width=5,height=5)
image(x=$cnvec,y=$cevec,z=$mpersist_surf,col=pal,xlab='cn',ylab='ce')
contour(x=$cnvec,y=$cevec,z=$mpersist_surf, add = TRUE, drawlabels = TRUE)
dev.off()
"""


#x ~ row of z
#y ~ columns of z
filename = "figures/extvpers_cncecp.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(9,'YlGnBu')
pdf($namespace,width=5,height=5)
plot($(vec(mextrate)),$(vec(mpersistance)),pch='.',col='#00000010')
dev.off()
"""

