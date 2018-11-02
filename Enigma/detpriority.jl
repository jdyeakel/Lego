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
keepstates = Array{Int64}(undef,0);
for i=1:nstates
    speciesobjects = states[i];
    adjacencymatrix = a_b[speciesobjects,speciesobjects];
    g = DiGraph(adjacencymatrix');
    paths = gdistances(g,1)
    if maximum(paths) < N+1;
        push!(keepstates,i);
    end
    
    #NOTE: take out states without complete set of species/object pairs
    
end
possiblestates = states[keepstates];

colonizers = potcol(sp_v,int_id,speciesobjects,a_b,n_b0,0,1);
