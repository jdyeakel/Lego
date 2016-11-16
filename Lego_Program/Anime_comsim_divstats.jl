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

@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")


##############################
# pr(need) vs. need_threshold
##############################


#Read-only variables
#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0;

tmax = 500;
reps=20;

needvec = collect(0.5:0.5:5.0);
n_thresh_vec = collect(0.0:0.05:0.5);

ln = length(needvec);
ltn = length(n_thresh_vec);

SPRICH = Array(Array{Int64},ln,ltn);
RICH =  Array(Array{Int64},ln,ltn);

num_play = 500;
ppweight = 1/4;
sim=false;
par=true;

for i=1:ln
  for j=1:ltn
    
    n_thresh = n_thresh_vec[j]
    trophicload = 2;
    
    #Establish community template
    probs = [
    p_n=needvec[i]/num_play,
    p_a=0.01,
    p_m= 1/num_play, #needvec[i]/num_play, #1/num_play,
    p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
    ]
    
    int_m, sprich, rich, conn, comgen, ext_prim, ext_sec = repsim(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight,sim,par);
    
    SPRICH[i,j] = sprich;
    RICH[i,j] = rich;
    
    
  end #end ltn
end #end ln

save("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prn_nt.jld","SPRICH",SPRICH,"RICH",RICH);


#Load library
d = load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prn_nt.jld");
SPRICH = d["SPRICH"];
RICH = d["RICH"];
  
richss = Array{Float64}(ln,ltn);
richsd = Array{Float64}(ln,ltn);
richcv = Array{Float64}(ln,ltn);
for i=1:ln
  for j=1:ltn
    reprichss = Array{Float64}(reps);
    reprichsd = Array{Float64}(reps);
    reprichcv = Array{Float64}(reps);
    for r=1:reps
      traj = copy(SPRICH[i,j][tmax-300:tmax,r])
      reprichss[r] = mean(traj);
      reprichsd[r] = std(traj);
      reprichcv[r] = std(traj)/mean(traj);
    end
    richss[i,j] = mean(reprichss);
    richsd[i,j] = mean(reprichsd);
    richcv[i,j] = mean(reprichcv);
  end
end


#Plotting the steady state mean
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/"
R"""
#pdf(paste($namespace,'fig_sen_prn_nt_mean.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral')); 
filled.contour(x=($needvec/$num_play),y=$n_thresh_vec,z=$richss,xlab='Pr(need)',ylab='n_thresh',key.title=title(main = 'mean'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
#dev.off()
"""

#Plotting steady state variance (SD)
R"""
#pdf(paste($namespace,'fig_sen_prn_nt_sd.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral')); 
filled.contour(x=($needvec/$num_play),y=$n_thresh_vec,z=$richsd,xlab='Pr(need)',ylab='n_thresh',key.title=title(main = 'sd'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
#dev.off()
"""

#Plotting steady state CV
R"""
#pdf(paste($namespace,'fig_sen_prn_nt_cv.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral')); 
filled.contour(x=($needvec/$num_play),y=$n_thresh_vec,z=$richcv,xlab='Pr(need)',ylab='n_thresh',key.title=title(main = 'cv'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
#dev.off()
"""



##############################
# pr(make) vs. need_threshold
##############################

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

@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")




#Read-only variables
#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0;

tmax = 500;
reps=20;

makevec = collect(0.5:0.5:5.0);
n_thresh_vec = collect(0.0:0.05:0.5);

ln = length(makevec);
ltn = length(n_thresh_vec);

SPRICH = Array(Array{Int64},ln,ltn);
RICH =  Array(Array{Int64},ln,ltn);

num_play = 500;
ppweight = 1/4;
sim=true;
par=true;

for i=1:ln
  for j=1:ltn
    
    n_thresh = n_thresh_vec[j]

    #Establish community template
    probs = [
    p_n=0.008, #This is chosen because there is rich behavior here over n_thresh
    p_a=0.01,
    p_m= makevec[i]/num_play, #needvec[i]/num_play, #1/num_play,
    p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
    ]
    
    int_m, sprich, rich, conn, comgen, ext_prim, ext_sec = repsim(num_play,reps,tmax,a_thresh,n_thresh,rate_col,probs,ppweight,sim,par);
    
    SPRICH[i,j] = sprich;
    RICH[i,j] = rich;
    
    
  end #end ltn
end #end ln

save("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prm_nt.jld","SPRICH",SPRICH,"RICH",RICH);


#Load library
d = load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prm_nt.jld");
SPRICH = d["SPRICH"];
RICH = d["RICH"];
  
richss = Array{Float64}(ln,ltn);
richsd = Array{Float64}(ln,ltn);
richcv = Array{Float64}(ln,ltn);
for i=1:ln
  for j=1:ltn
    reprichss = Array{Float64}(reps);
    reprichsd = Array{Float64}(reps);
    reprichcv = Array{Float64}(reps);
    for r=1:reps
      traj = copy(SPRICH[i,j][tmax-300:tmax,r])
      reprichss[r] = mean(traj);
      reprichsd[r] = std(traj);
      reprichcv[r] = std(traj)/mean(traj);
    end
    richss[i,j] = mean(reprichss);
    richsd[i,j] = mean(reprichsd);
    richcv[i,j] = mean(reprichcv);
  end
end


#Plotting the steady state mean
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/"
R"""
pdf(paste($namespace,'fig_sen_prm_nt_mean.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral')); 
filled.contour(x=($needvec/$num_play),y=$n_thresh_vec,z=$richss,xlab='Pr(make)',ylab='n_thresh',key.title=title(main = 'mean'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
dev.off()
"""

#Plotting steady state variance (SD)
R"""
pdf(paste($namespace,'fig_sen_prm_nt_sd.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral')); 
filled.contour(x=($needvec/$num_play),y=$n_thresh_vec,z=$richsd,xlab='Pr(make)',ylab='n_thresh',key.title=title(main = 'sd'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
dev.off()
"""

#Plotting steady state CV
R"""
pdf(paste($namespace,'fig_sen_prm_nt_cv.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral')); 
filled.contour(x=($needvec/$num_play),y=$n_thresh_vec,z=$richcv,xlab='Pr(make)',ylab='n_thresh',key.title=title(main = 'cv'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
dev.off()
"""
