using Distributed
using DataFrames
using Images

using AxisArrays
using Combinatorics
using LinearAlgebra
# @everywhere using Distributed
using SharedArrays
using SparseArrays
# using DataFrames
using Distributions
# using Images
using SpecialFunctions
using LightGraphs
using RCall
# @everywhere using HDF5
using JLD2

@everywhere using AxisArrays
@everywhere using Combinatorics
@everywhere using LinearAlgebra
# @everywhere using Distributed
@everywhere using SharedArrays
@everywhere using SparseArrays
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using Images
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
# @everywhere using HDF5
@everywhere using JLD2

if homedir() == "/home/z840"

    @everywhere include("$(homedir())/2014_Lego/Enigma/src/smartpath.jl")

    #Interaction matrix
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/intmatrixv3.jl")

		#Alternative 4x4 Interaction matrix
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/intmatrixv4.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/intmatrixv4_unique.jl")

    #Community dynamics
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/preamble_defs.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/assembly.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/assembly_delayedobjects.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/assemblystate.jl")

    #Analysis Calculations
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/structure.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/dynstructure.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/sortassembly.jl")

    #Analysis functions
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/trophicalc2.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/roverlap.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/potcol.jl")

    @everywhere include("$(homedir())/2014_Lego/Enigma/src/nichemodelweb.jl")
    @everywhere include("$(homedir())/2014_Lego/Enigma/src/proposal.jl")


else

    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/smartpath.jl")

    #Interaction matrix
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/intmatrixv3.jl")

		#Alternative 4x4 Interaction matrix
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/intmatrixv4.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/intmatrixv4_unique.jl")

    #Community dynamics
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/preamble_defs.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly_delayedobjects.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assemblystate.jl")


    #Analysis Calculations
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/structure.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/dynstructure.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/sortassembly.jl")

    #Analysis functions
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/trophicalc2.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/roverlap.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/potcol.jl")

    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/nichemodelweb.jl")
    @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/proposal.jl")
end
