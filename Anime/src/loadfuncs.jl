@everywhere using Distributions
#@everywhere using Gadfly
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD

#Interaction matrix
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/build_template_species.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/intmatrix.jl")

#Community dynamics
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/colext.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/preamble_defs.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/assembly.jl")

#Analysis Calculations
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/structure.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/dynstructure.jl")

#Analysis functions
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/trophicalc2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/pairbin.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/trophicwidth.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/Jmatrix.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/PSWebs.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/movingaverage.jl")
