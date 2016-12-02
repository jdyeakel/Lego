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
reps=50;

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
POT_COL = Array(Array{Int64},lm);


for i=3:3
  println("i=",i)
  #Establish community template
  probs = [
  p_n=0.004,
  p_a=0.01,
  p_m= makevec[i], #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]

  int_mv, sprich, rich, conn, comgen, ext_prim, ext_sec, pot_col = repsimint(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  SPRICH[i] = sprich;
  RICH[i] = rich;
  EXTINCTIONS_PRIM[i] = ext_prim;
  EXTINCTIONS_SEC[i] = ext_sec;
  POT_COL[i] = pot_col;
end

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_make/rich_prm_tl",trophicload,".jld");
save(namespace,"makevec",makevec,"SPRICH",SPRICH,"RICH",RICH,"EXTINCTIONS_PRIM",EXTINCTIONS_PRIM,"EXTINCTIONS_SEC",EXTINCTIONS_SEC,"POT_COL",POT_COL);

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
#makevec = d["makevec"];
SPRICH = d["SPRICH"];
RICH = d["RICH"];
EXTINCTIONS_PRIM = d["EXTINCTIONS_PRIM"];
EXTINCTIONS_SEC = d["EXTINCTIONS_SEC"];
POT_COL = d["POT_COL"];

startpt = 500;
n=1000;
lm = length(SPRICH);
reps = length(SPRICH[1][1,:]);
r2rich = Array{Float64}(lm,reps);
r2ext = Array{Float64}(lm,reps);
corrich = Array{Float64}(lm,reps);
corext = Array{Float64}(lm,reps);
cordiff = Array{Float64}(lm,reps);
cordiffbin = Array{Float64}(lm,reps);
corcol = Array{Float64}(lm,reps);
makevec = collect(0:0.005:0.03);

#Analysis
for i=1:lm
  sprich_i = SPRICH[i];
  rich_i = RICH[i];
  extprim_i = EXTINCTIONS_PRIM[i];
  extsec_i = EXTINCTIONS_SEC[i];
  potcol_i = POT_COL[i];
  for r=1:reps
    sprich=sprich_i[:,r];
    rich = rich_i[:,r];
    obrich = rich.-sprich;
    extprim = extprim_i[:,r];
    extsec = extsec_i[:,r];
    potcol_traj = potcol_i[:,r];
    tmax = length(sprich);
    #Grab random time points across the trajectory
    #Such that rdsample = t, and rdsample+1 = t+1
    wipeout = find(x->x==0,rich);
    timesample = collect(startpt:tmax-1);
    ts_exist = zeros(Int64,0);
    for j=1:length(timesample)
      if in(timesample[j],wipeout)==false
        append!(ts_exist,timesample[j]);
      end
    end
    rdsample = sample(ts_exist,n,replace=false);
    #Objects per Total Richness at time t
    obj_ratio = obrich[rdsample] #./ rich[rdsample];
    #Species per Total Richness at time t+1
    sprich_r = sprich[rdsample+1] ./ rich[rdsample+1];
    #Number of extinctions at time t+1 per Species at time t
    ext_ratio = (extprim[rdsample+1]+extsec[rdsample+1]) #./ rich[rdsample];
    #Change in species from t to t+1 per Total richness at time t
    spdiff = (sprich[rdsample+1]-sprich[rdsample]) #./ rich[rdsample];
    #Potential colonizers at timestep t
    potcol = potcol_traj[rdsample+1];
    
    nonext = find(x->x!=0,ext_ratio);
    nondiff = find(x->x!=0,spdiff);
    #nz = length(nonzero);
    #R"lmrich=lm($(sprich_r[nonzero]) ~ $(obj_ratio[nonzero]));
    #r2r = summary(lmrich)$adj.r.squared";
    #r2rich[i,r] = @rget r2r;
    #corrich[i,r] = cor(obj_ratio[nonzero],sprich_r[nonzero]);
    corext[i,r] = cor(obj_ratio[nonext],ext_ratio[nonext]);
    cordiff[i,r] = cor(obj_ratio[nondiff],spdiff[nondiff]);
    corcol[i,r] = cor(obj_ratio,potcol);
    R"lmext=lm($(ext_ratio[nonext]) ~ $(obj_ratio[nonext]));
    r2e = summary(lmext)$adj.r.squared";
    R"lmdiff=lm($(spdiff[nondiff]) ~ $(obj_ratio[nondiff]));
    r2e = summary(lmdiff)$adj.r.squared";
    R"lmcol = lm($(potcol) ~ $(obj_ratio))";
    r2ext[i,r] = @rget r2e;
  end
end
#

#boxplot(t($corrich),names=$makevec,boxwex=0.5,xlab='pr(m)',ylab='Corr Ob(t)/R(t) vs Sp(t+1)/R(t+1)',col=cols[2])

#Correlation
namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_corrobext_tl",trophicload,".pdf");
R"""
library(RColorBrewer)
cols <- brewer.pal(3,'Set1')
pdf($namespace,height=6,width=12)
par(mfrow=c(1,2))
boxplot(t($corext),names=$makevec,boxwex=0.5,xlab='pr(m)',ylab='Corr Ob(t) vs Ext(t+1)',col=cols[2])
lines(seq(0,10,length.out=5),seq(0,0,length.out=5),lty=3)
boxplot(t($corcol),names=$makevec,boxwex=0.5,xlab='pr(m)',ylab='Corr Ob(t) vs Col(t+1)',col=cols[2])
lines(seq(0,10,length.out=5),seq(0,0,length.out=5),lty=3)
dev.off()
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
plot($(obj_ratio[nonext]),$(ext_ratio[nonext]))
abline(lmext)
plot($(obj_ratio),$(potcol))
abline(lmcol)
"""

R"plot($sprich,type='l',ylim=c(0,max($rich)));
lines($rich,lty=3)"
