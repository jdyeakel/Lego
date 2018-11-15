#RUN REMOTELY
if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");

# RUN LOCALLY
using Distributed
@everywhere using Combinatorics
@everywhere using LinearAlgebra
# @everywhere using Distributed
@everywhere using SharedArrays
@everywhere using SparseArrays
@everywhere using Distributions
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
# @everywhere using HDF5
@everywhere using JLD
#Interaction matrix
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/intmatrixv3.jl")
#Community dynamics
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/preamble_defs.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly.jl")
#Analysis Calculations
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/structure.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/dynstructure.jl")
#Analysis functions
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/trophicalc2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/roverlap.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/potcol.jl")


S = 5;

# S = 400;
probs = (
p_n=0.1,
p_a=0.3
# p_n = 0.02,
# p_a = 0.02
);
#expected objects per species
lambda = 0.0;

MaxN = convert(Int64,floor(S + S*lambda));

int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

a_b,
n_b,
i_b,
m_b,
n_b0,
sp_v,
int_id = preamble_defs(int_m);

N = size(int_m)[1];

#Define all possible states
states = collect(combinations(collect(1:N)));
states = states[N+1:length(states)];
states = states[in.(1,states)];
#how many states
nstates = length(states);

#Determine impossible states
passtest = zeros(Int64,nstates) .+ 1;
keepstates = Array{Int64}(undef,0);
for i=1:nstates
    speciesobjects = states[i];
    adjacencymatrix = a_b[speciesobjects,speciesobjects] .+ m_b[speciesobjects,speciesobjects]';
    g = DiGraph(adjacencymatrix');
    paths = gdistances(g,1)
    if maximum(paths) < N+1;
        passtest[i] *= 1;
    else
        passtest[i] *= 0;
    end
    
    #NOTE: take out states without complete set of species/object pairs
    observedspecies = speciesobjects[speciesobjects .<= S];
    expectedobjects = findall(!iszero,vec(sum(m_b[observedspecies,:],dims=1)));
    
    observedobjects = setdiff(speciesobjects,observedspecies);
    expectedspecies = findall(!iszero,vec(sum(m_b[:,observedobjects],dims=2)));
    
    if observedspecies == expectedspecies || observedobjects == expectedobjects
        passtest[i] *= 1;
    else
        passtest[i] *= 0;
    end
    
    
end

possiblestates = states[findall(!iszero,passtest)];
statespace = sum(passtest)/nstates
lstate = length(possiblestates);

# SparseArray
transm = spzeros(lstate,lstate);
for i=1:lstate
    # print(string(i,'_'))
    statei = copy(possiblestates[i]);
    deleteat!(statei,1)
    colonizers = potcol(sp_v,int_id,statei,a_b,n_b0,0,1);
    # newstates = Array{Array}(undef,length(colonizers));
    newstatesloc = Array{Int64}(undef,length(colonizers));
    for j=1:length(colonizers)
        newstates = sort([1;statei;colonizers[j]]);
        newstatesloc[j] = findall(x->x==newstates,possiblestates)[1];
    end
    transm[i,newstatesloc] .= 1.0;
end

transg = DiGraph(transm);

R"""
library(igraph)
library(RColorBrewer)
#pdf($namespace,width=6,height=5)
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(Array(transm)));
coords <- layout_(fw_g, as_tree(circular=FALSE))
plot(fw_g,layout=coords,vertex.size=2,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA,vertex.frame.color=NA) 
#dev.off()
"""

R"""
g = graph.adjacency($(Array(transm)));
plot(g,layout=layout_as_tree(g,circular=TRUE))
"""

,vertex.size=2,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA,vertex.frame.color=NA