using Distributed

@everywhere using LinearAlgebra
@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
# @everywhere using HDF5
@everywhere using JLD

#Interaction matrix
@everywhere include("$(homedir())/2014_Lego/Enigma/src/intmatrixv3.jl")

#Community dynamics
@everywhere include("$(homedir())/2014_Lego/Enigma/src/preamble_defs.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/assembly.jl")

#Analysis Calculations
@everywhere include("$(homedir())/2014_Lego/Enigma/src/structure.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/dynstructure.jl")

#Analysis functions
@everywhere include("$(homedir())/2014_Lego/Enigma/src/trophicalc2.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/roverlap.jl")
@everywhere include("$(homedir())/2014_Lego/Enigma/src/potcol.jl")
