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


species = Array{Float64}(undef,annealtime,2);
conn = Array{Float64}(undef,annealtime,2);
mdegree = Array{Float64}(undef,annealtime,2);
stdindegree = Array{Float64}(undef,annealtime,2);
stdoutdegree = Array{Float64}(undef,annealtime,2);
#starting error
# global errvec_old = [5];
#starting temperature

mtemp = Array{Float64}(undef,annealtime);
tempvec = Array{Float64}(undef,annealtime);
error = Array{Float64}(undef,annealtime);
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

#Search for parameters that match niche model
annealtime = 50;

#The range of values to explore
cnrange = collect(0.1:(10-.1)/999:10);
cerange = collect(0.1:(10-.1)/999:10);
cprange = collect(0.1:(10-.1)/999:10);
p_nrange = collect(0.0005:(0.003-0.0005)/999:0.003);
p_arange = collect(0.005:(0.03-0.005)/999:0.03);


# cn = cnvec[208];
# ce = cevec[134];
# cp = cpvec[92];
# p_n = p_nvec[11];
# p_a = p_avec[75];
# errmean_old = 5.;
# temperature = 1.0;
# tic = 1;

cn_init = 3.14;
ce_init = 1.41;
cp_init = 1.0;
p_n_init = 0.002;
p_a_init = 0.01;

let cn = cnrange[findall(x->x==minimum(abs.(cnrange .- cn_init)),abs.(cnrange .- cn_init))[1]], 
    ce = cerange[findall(x->x==minimum(abs.(cerange .- ce_init)),abs.(cerange .- ce_init))[1]], 
    cp = cprange[findall(x->x==minimum(abs.(cprange .- cp_init)),abs.(cprange .- cp_init))[1]],
    p_n = p_nrange[findall(x->x==minimum(abs.(p_nrange .- p_n_init)),abs.(p_nrange .- p_n_init))[1]],
    p_a = p_arange[findall(x->x==minimum(abs.(p_arange .- p_a_init)),abs.(p_arange .- p_a_init))[1]],
    errmean_old = 1.,
    temperature = 1.0,
    tic = 1

    # global cn = 4.195216435170081;
    # global ce = 1.8682437375935583;
    # global cp = 0.6511860914279262;
    # global p_n = 0.0030274490082625636;
    # global p_a =  0.014690692329311795;


    for r=1:annealtime

        tempvec[r] = copy(temperature);
        
        #Select proposed value from predefined range with uniform selection over a smaller range determined by temperature
        cn_prop = proposal(temperature/2,cn,cnrange);
        ce_prop = proposal(temperature/2,ce,cerange);
        cp_prop = proposal(temperature/2,cp,cprange);
        p_n_prop = proposal(temperature/2,p_n,p_nrange);
        p_a_prop = proposal(temperature/2,p_a,p_arange);
        
        
        S = 200;
        maxits = 2000;
        SOprobs = (
        p_n=p_n_prop,
        p_a=p_a_prop
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
                athresh,nthresh,maxits,cn_prop,ce_prop,cp_prop);

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
        mdegree[r,2] = mean(jmdegree[vec(findall(!isnan,jmdegree))]);
        stdindegree[r,2] = mean(jstdindegree[vec(findall(!isnan,jstdindegree))]);
        stdoutdegree[r,2] = mean(jstdoutdegree[vec(findall(!isnan,jstdoutdegree))]);


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
        # err_sp = sqrt((species[r,2] - species[r,1])^2)/mean(species[r,1]); #/std(ispecies);
        # err_conn = sqrt((conn[r,2] - conn[r,1])^2)/mean(conn[r,1]); #/std(iconn);
        # err_md = sqrt((mdegree[r,2] -mdegree[r,1])^2)/mean(mdegree[r,1]); #/std(imdegree);
        # err_sdin = sqrt((stdindegree[r,2] - stdindegree[r,1])^2)/mean(stdindegree[r,1]); #/std(istdindegree);
        # err_sdout = sqrt((stdoutdegree[r,2] - stdoutdegree[r,1])^2)/mean(stdoutdegree[r,1]); #/std(istdoutdegree);

        errvec = sqrt((mdegree[r,2] -mdegree[r,1])^2); #/std(imdegree);
        # errvec = [err_sp,err_conn,err_md,err_sdin,err_sdout];
        
        # errvec = [err_md,err_sdin,err_sdout]

        error[r] = copy(errvec);
        errmean = mean(errvec);
        
        

        errscore[r] = errmean;

        prob_accept = exp(-(errmean-errmean_old)/temperature);
        rdraw = rand();
        
        if errmean < errmean_old
            #Accept and lower the temperature
            cn = copy(cn_prop);
            ce = copy(ce_prop);
            cp = copy(cp_prop);
            p_n = copy(p_n_prop);
            p_a = copy(p_a_prop);
            
            temperature = temperature*(errmean/errmean_old);
        end
        # elseif rdraw < prob_accept
        #     #Accept and lower the temperature
        #     cn = copy(cn_prop);
        #     ce = copy(ce_prop);
        #     cp = copy(cp_prop);
        #     p_n = copy(p_n_prop);
        #     p_a = copy(p_a_prop);
        # 
        #     temperature = temperature*(errmean/errmean_old);
        # end
        
        println(string(r,": error=",round(errmean,digits=3),"; temp=",round(temperature,digits=3),"; P=",round(prob_accept,digits=3)))
        
        cnvec[r] = cn;
        cevec[r] = ce;
        cpvec[r] = cp;
        pnvec[r] = p_n;
        pevec[r] = p_a;
        
        #update error
        errvec_old = copy(errvec);
        errmean_old = copy(errmean);


    end
end

filename = "figures/niche/annealingtemp.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=5,width=6)
plot($(tempvec),type='l',ylim=c(10^(-4),0.5),log='y',col=pal[1])
"""
R"dev.off()"

filename = "figures/niche/errorscore.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=5,width=6)
plot($(error),type='l',col=pal[1])
"""
R"dev.off()"


filename = "figures/niche/annealingerror.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=5,width=6)
plot($(error[:,1]),type='l',ylim=c(10^(-1),1),log='y',col=pal[1])
"""
for i=2:5
    R"""
    lines($(error[:,i]),col=pal[$i])
    """
end
R"lines($errscore,col='black')"
R"dev.off()"



filename = "figures/niche/probaccept.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=5,width=6)
plot($(tempvec[2:size(tempvec)[1]]),$(prob_accept),ylim=c(0,2))
"""
R"dev.off()"



filename = "figures/niche/vec.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=5,width=6)
plot($(cnvec),type='l',ylim=c(0,3),lty=1)
lines($cevec,lty=2)
lines($cpvec,lty=3)
"""
R"dev.off()"
