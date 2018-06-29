@everywhere using Distributions
@everywhere using SpecialFunctions
@everywhere using LightGraphs
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD
# 
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/intmatrixv3.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/preamble_defs.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly.jl")


S = 400;

tmax = 10;

# S = 400;
probs = [
p_n=0.003,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
lambda = 0.5;
athresh = 0;
nthresh = 0.5;
MaxN = convert(Int64,floor(S + S*lambda));

sprich = Array{Int64}(tmax);
int_m = Array{Char}();
tp_m = Array{Int64}();
tind_m = Array{Int64}();
prim_ext = Array{Int64}(tmax);
sec_ext = Array{Int64}(tmax);
status = Array{Int64}(tmax);
lpot_col = Array{Int64}(tmax);
CID = Array{Bool}(tmax,MaxN)*false;

int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

a_b,
n_b,
i_b,
m_b,
n_b0,
sp_v,
int_id = preamble_defs(int_m);

@time sprich,rich,clock = assembly(
    int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
    athresh,nthresh,tmax);

R"plot($sprich,type='l')"
R"lines($(rich .- sprich),col='gray')"


