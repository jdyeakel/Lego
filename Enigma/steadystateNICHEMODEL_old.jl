if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

filename = "data/steadystate/sim_settings.jld";
namespace = smartpath(filename);
@load namespace reps S maxits athresh nthresh lambda SSprobs SOprobs OOprobs;

species = SharedArray{Float64}(reps,2);
conn = SharedArray{Float64}(reps,2);
stdindegree = SharedArray{Float64}(reps,2);
stdoutdegree = SharedArray{Float64}(reps,2);

mpath = SharedArray{Float64}(reps,2);
stdpath = SharedArray{Float64}(reps,2);
maxpath = SharedArray{Float64}(reps,2);


@sync @distributed for r=1:reps
    #Read in the interaction matrix


    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;


    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;


    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);


    Aenigma = a_b[CID[:,maxits],CID[:,maxits]];
    Senigma = size(Aenigma)[1];
    Aenigma = Aenigma[2:Senigma,2:Senigma];
    Senigma -= 1;


    Cenigma = sum(Aenigma)/(Senigma^2);

    #NICHE MODEL COMPARISON
    Aniche, n = nichemodelweb(Senigma,Cenigma);

    Aniche = Array{Int64}(Aniche);

    Sniche = size(Aniche)[1];
    Cniche = sum(Aniche)/Sniche^2;

    species[r,:] = [Sniche,Senigma];
    #Connectance
    conn[r,:] = [Cniche,Cenigma];

    #indegree
    stdindegree[r,1] = std(sum(Aniche,dims=2));
    stdindegree[r,2] = std(sum(Aenigma,dims=2));

    #outdegree
    stdoutdegree[r,1] = std(sum(Aniche,dims=1));
    stdoutdegree[r,2] = std(sum(Aenigma,dims=1));

    #mean chainlength
    gniche = DiGraph(Aniche');
    genigma = DiGraph(Aenigma');
    paths_niche = Array{Int64}(undef,Sniche,Sniche);
    paths_enigma = Array{Int64}(undef,Senigma,Senigma);
    for i=1:maximum([Sniche,Senigma])
        if i <= Sniche
            paths_niche[:,i] = gdistances(gniche,i);
        end
        if i <= Senigma
            paths_enigma[:,i] = gdistances(genigma,i);
        end
    end
    paths_nichepos = findall(x->x<Sniche^2,paths_niche);
    paths_enigmapos = findall(x->x<Senigma^2,paths_enigma);

    mpath[r,1] = mean(paths_niche[paths_nichepos]);
    mpath[r,2] = mean(paths_enigma[paths_enigmapos]);

    stdpath[r,1] = std(paths_niche[paths_nichepos]);
    stdpath[r,2] = std(paths_enigma[paths_enigmapos]);

    maxpath[r,1] = maximum(paths_niche[paths_nichepos]);
    maxpath[r,2] = maximum(paths_enigma[paths_enigmapos]);

    #omnivory : fraction of species that consume >= 2 things with diff chain lengths

    #SD chainlength
end

filename = "figures/comparison.pdf"
namespace = smartpath(filename);
R"""
pdf($namespace,height=5,width=6)
boxplot($stdindegree)
dev.off()
"""

filename = "figures/foodweb.pdf"
namespace = smartpath(filename);
adjmatrix = Aniche;
R"""
par(mfrow=c(1,2))
library(igraph)
library(RColorBrewer)
pdf($namespace,width=6,height=5)
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(adjmatrix'));
plot(fw_g,layout=layout_with_fr())
"""
adjmatrix = Aenigma;
R"""
fw_g <- graph.adjacency($(adjmatrix'));
plot(fw_g,layout=layout_with_fr())
dev.off()
"""


cid = findall(isodd,CID[:,maxits]);
deg,troph = structure(S,cid,sp_v,tind_m);
filename = "figures/foodweb.pdf";
namespace = smartpath(filename);
adjmatrix = Aenigma;
R"""
library(igraph)
library(RColorBrewer)
pdf($namespace,width=6,height=5)
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
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)))
#main=ecount(fw_g)/$(size(adjmatrix)[1])^2,
fw_ind <- graph.adjacency($(indmatrix[keepnodes,keepnodes]'));
#plot(fw_ind,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='red',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)),add=TRUE)
dev.off()
"""



if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

#Search for parameters that match niche model
annealtime = 50;
global cn = pi;
global ce = sqrt(2);
global cp = 1;
global p_n = 0.002;
global p_a = 0.02;
global tic = 1;

species = Array{Float64}(undef,annealtime,2);
conn = Array{Float64}(undef,annealtime,2);
mdegree = Array{Float64}(undef,annealtime,2);
stdindegree = Array{Float64}(undef,annealtime,2);
stdoutdegree = Array{Float64}(undef,annealtime,2);
global errmean_old = 100.;
global errvec_old = [100,100,100,100,100];
global temperature = [0.5,0.5,0.5,0.5,0.5];

mtemp = Array{Float64}(undef,annealtime);
tempvec = Array{Float64}(undef,annealtime,5);
error = Array{Float64}(undef,annealtime,5);
errscore = Array{Float64}(undef,annealtime);
cnvec = Array{Float64}(undef,annealtime);
cevec = Array{Float64}(undef,annealtime);
cpvec = Array{Float64}(undef,annealtime);
pnvec = Array{Float64}(undef,annealtime);
pevec = Array{Float64}(undef,annealtime);


# temperature = tempvec[50,:];
# global cn = cnvec[50];
# global ce = cevec[5];
# global cp = cpvec[50];
# global p_n = pnvec[50];
# global p_a = pevec[50];
# global tic = 1;
#


for r=1:annealtime

    tempvec[r,:] = temperature;


    altdist = Normal.(0,temperature);
    altvec = rand.(altdist);

    global cn = cn*(1 + altvec[1]);
    global ce = ce*(1 + altvec[2]);
    global cp = cp*(1 + altvec[3]);

    global p_n=p_n*(1 + altvec[4]);
    global p_a=p_a*(1 + altvec[5]);

    S = 200;
    maxits = 2000;
    SOprobs = (
    p_n=p_n,
    p_a=p_a
    );
    SSmult = 1.0; OOmult = 0.0;
    SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
    OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


    #expected objects per species
    lambda = 0;
    athresh = 0;
    nthresh = 1.0;
    MaxN = convert(Int64,floor(S + S*lambda));

    enigmareps=50;
    jspecies = SharedArray{Float64}(enigmareps);
    jconn = SharedArray{Float64}(enigmareps);
    jmdegree = SharedArray{Float64}(enigmareps);
    jstdindegree = SharedArray{Float64}(enigmareps);
    jstdoutdegree = SharedArray{Float64}(enigmareps);
    @sync @distributed for j = 1:enigmareps
        int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

        a_b,
        n_b,
        i_b,
        m_b,
        n_b0,
        sp_v,
        int_id = preamble_defs(int_m);

        sprich,rich,clock,CID = assembly(
            int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
            athresh,nthresh,maxits,cn,ce,cp);

        Aenigma = a_b[CID[:,maxits],CID[:,maxits]];
        Senigma = size(Aenigma)[1];
        Aenigma = Aenigma[2:Senigma,2:Senigma];
        Senigma -= 1;
        Cenigma = sum(Aenigma)/(Senigma^2);

        jspecies[j] = Senigma;
        jconn[j] = Cenigma;
        jmdegree[j] = mean(sum(Aenigma,dims=2));
        jstdindegree[j] = std(sum(Aenigma,dims=2));
        jstdoutdegree[j] = std(sum(Aenigma,dims=1));
    end



    species[r,2] = mean(jspecies[vec(findall(!isnan,jspecies))]);
    conn[r,2] = mean(jconn[vec(findall(!isnan,jconn))]);
    mdegree[r,2] = mean(jmdegree);
    stdindegree[r,2] = mean(jstdindegree);
    stdoutdegree[r,2] = mean(jstdoutdegree);


    #Simulate a bunch of nichemodelwebs



    nichereps = 100;
    ispecies = Array{Float64}(undef,nichereps);
    iconn = Array{Float64}(undef,nichereps);
    imdegree = Array{Float64}(undef,nichereps);
    istdindegree = Array{Float64}(undef,nichereps);
    istdoutdegree = Array{Float64}(undef,nichereps);

    for i=1:nichereps
        #NICHE MODEL COMPARISON
        #NOTE: C CANNOT BE EQUAL TO OR GREATER THAN 0.5
        #NOTE Violates beta parameters!
        Aniche, n = nichemodelweb(Int64(floor(species[r,2])),conn[r,2]);

        #make measurements
        Aniche = Array{Int64}(Aniche);

        Sniche = size(Aniche)[1];
        Cniche = sum(Aniche)/Sniche^2;

        ispecies[i] = Sniche;
        #Connectance
        iconn[i] = Cniche;

        #mean degree
        imdegree[i] = mean(sum(Aniche,dims=2));

        #indegree
        istdindegree[i] = std(sum(Aniche,dims=2));

        #outdegree
        istdoutdegree[i] = std(sum(Aniche,dims=1));
    end

    #ignore NAN
    
    species[r,1] = mean(ispecies);

    #Connectance
    conn[r,1] = mean(iconn);

    #mean degree
    mdegree[r,1] = mean(imdegree);


    #indegree
    stdindegree[r,1] = mean(istdindegree);


    #outdegree
    stdoutdegree[r,1] = mean(istdoutdegree);


    #Calculate error
    z_sp = sqrt((species[r,2] - mean(ispecies))^2); #/std(ispecies);
    z_conn = sqrt((conn[r,2] - mean(iconn))^2); #/std(iconn);
    z_md = sqrt((mdegree[r,2] - mean(imdegree))^2); #/std(imdegree);
    z_sdin = sqrt((stdindegree[r,2] - mean(istdindegree))^2); #/std(istdindegree);
    z_sdout = sqrt((stdoutdegree[r,2] - mean(istdoutdegree))^2); #/std(istdoutdegree);

    global zvec = [z_sp,z_conn,z_md,z_sdin,z_sdout];

    zmean = mean(zvec);

    #temperature goes down as zvec gets smaller
    # global temperature = temperature .* (zvec ./ zvec_old);

    #Only lower temperature if
    for i=1:length(zvec)
        if zvec[i] < zvec_old[i]
            temperature[i] = temperature[i] * (zvec[i]/zvec_old[i])
        end
    end

    mtemp[r] = mean(temperature);


    global zmean_old = copy(zmean);
    global zvec_old = copy(zvec);

    cnvec[r] = cn;
    cevec[r] = ce;
    cpvec[r] = cp;
    pnvec[r] = SOprobs.p_n;
    pevec[r] = SOprobs.p_a;

    zscore[r] = abs(zmean);


    println(string(r,": errscore=",errscore[r]))
end

filename = "figures/niche/annealingtempnew.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=5,width=6)
plot($(tempvec[:,1]),type='l',ylim=c(0.00000000001,0.5),log='y',col=pal[1])
"""
for i=2:5
    R"""
    lines($(tempvec[:,i]),col=pal[$i])
    """
end
R"dev.off()"
