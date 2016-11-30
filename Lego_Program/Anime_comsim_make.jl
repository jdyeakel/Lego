using Distributions
#using Gadfly
using RCall
using HDF5
using JLD

@everywhere using Distributions
#@everywhere using Gadfly
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsim.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsimint.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")




#Read-only variables
#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0;
n_thresh = 0.2;

tmax = 5000;
reps=100;

num_play = 500;
ppweight = 1/4;
trophicload=2;

#Search over pr_m
makevec = collect(0:0.005:0.03);
lm = length(makevec);

SPRICH = Array(Array{Int64},lm);
RICH =  Array(Array{Int64},lm);
EXTINCTIONS_PRIM = Array(Array{Int64},lm);
EXTINCTIONS_SEC = Array(Array{Int64},lm);


for i=1:lm
  println("i=",i)
  #Establish community template
  probs = [
  p_n=0.004,
  p_a=0.01,
  p_m= makevec[i], #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]

  int_mv, sprich, rich, conn, comgen, ext_prim, ext_sec = repsimint(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  SPRICH[i] = sprich;
  RICH[i] = rich;
  EXTINCTIONS_PRIM[i] = ext_prim;
  EXTINCTIONS_SEC[i] = ext_sec;
end

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_make/rich_prm_tl",trophicload,".jld");
save(namespace,"SPRICH",SPRICH,"RICH",RICH,"EXTINCTIONS_PRIM",EXTINCTIONS_PRIM,"EXTINCTIONS_SEC",EXTINCTIONS_SEC);
quit();


#########################
#########################
using Distributions
#using Gadfly
using RCall
using HDF5
using JLD

#Load library
trophicload = 2;
namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_make/rich_prm_tl",trophicload,".jld");
d = load(namespace);
SPRICH = d["SPRICH"];
RICH = d["RICH"];
EXTINCTIONS_PRIM = d["EXTINCTIONS_PRIM"];
EXTINCTIONS_SEC = d["EXTINCTIONS_SEC"];


startpt = 300;
n=3000;
lm = length(SPRICH);
reps = length(SPRICH[1][1,:]);
r2rich = Array{Float64}(lm,reps);
r2ext = Array{Float64}(lm,reps);
corrich = Array{Float64}(lm,reps);
corext = Array{Float64}(lm,reps);
#Analysis
for i=1:lm
  sprich_i = SPRICH[i];
  rich_i = RICH[i];
  extprim_i = EXTINCTIONS_PRIM[i];
  extsec_i = EXTINCTIONS_SEC[i];
  for r=1:reps
    sprich=sprich_i[:,r];
    rich = rich_i[:,r];
    obrich = rich.-sprich;
    extprim = extprim_i[:,r];
    extsec = extsec_i[:,r];
    #find nonzero randomized elements
    tmax = length(sprich);
    rdsample = sample(collect(startpt:tmax-1),n,replace=false);
    obj_ratio = obrich[rdsample]./sprich[rdsample];
    sprich_r = sprich[rdsample+1]; #./ rich[rdsample+1];
    ext_ratio = (extprim[rdsample+1]+extsec[rdsample+1])./ sprich[rdsample];
    nonzero = find(x->x!=0,ext_ratio);
    nz = length(nonzero);
    R"lmrich=lm($(sprich_r[nonzero]) ~ $(obj_ratio[nonzero]));
    r2r = summary(lmrich)$adj.r.squared";
    r2rich[i,r] = @rget r2r;
    corrich[i,r] = cor(obj_ratio[nonzero],sprich_r[nonzero]);
    corext[i,r] = cor(obj_ratio[nonzero],ext_ratio[nonzero]);
    R"lmext=lm($(ext_ratio[nonzero]) ~ $(obj_ratio[nonzero]));
    r2e = summary(lmext)$adj.r.squared";
    r2ext[i,r] = @rget r2e;
  end
end

#Correlation
R"""
par(mfrow=c(1,2))
boxplot(t($corrich),names=$makevec,boxwex=0.5)
boxplot(t($corext),names=$makevec,boxwex=0.5)
"""

#RSquared
R"""
par(mfrow=c(1,2))
r2rich = $(r2rich);
boxplot(t(r2rich))
r2ext = $(r2ext);
boxplot(t(r2ext))
"""


#Plotting individual scenarios
R"""
par(mfrow=c(1,2))
plot($(obj_ratio[nonzero]),$(sprich_r[nonzero]))
abline(lmrich)
plot($(obj_ratio[nonzero]),$(ext_ratio[nonzero]))
abline(lmext)
"""

R"plot($sprich,type='l',ylim=c(0,max($rich)));
lines($rich,lty=3)"
