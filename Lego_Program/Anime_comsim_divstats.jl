# using Distributions
# #using Gadfly
# using RCall
# using HDF5
# using JLD
#
# @everywhere using Distributions
# #@everywhere using Gadfly
# @everywhere using RCall
# @everywhere using HDF5
# @everywhere using JLD
#
#
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsim.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsimintsingle.jl")
#
#
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
#
#
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonizesingle_func.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")
#
#
# ##############################
# # pr(need) vs. need_threshold
# ##############################
#
#
# #Read-only variables
# #Establish colonization and extinction rates
# rate_col = 1;
# #Establish thresholds
# a_thresh = 0;
#
# tmax = 2000;
# reps=25;
#
# needvec = collect(0.001:0.001:0.01);
# n_thresh_vec = collect(0.0:0.05:0.5);
#
# ln = length(needvec);
# ltn = length(n_thresh_vec);
#
# SPRICH = Array(Array{Int64},ln,ltn);
# RICH =  Array(Array{Int64},ln,ltn);
#
# S = 400;
# ppweight = 1/4;
# sim=false;
# par=true;
#
# for i=1:ln
#   for j=1:ltn
#     println("i=",i,"/",ln," j=",j,"/",ltn)
#
#     n_thresh = n_thresh_vec[j]
#     trophicload = 2;
#
#     #Establish community template
#     probs = [
#     p_n=needvec[i],
#     p_a=0.01,
#     p_m= 0.0008, #needvec[i]/num_play, #1/num_play,
#     p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
#     ]
#
#     sprich,
#     rich,
#     conn,
#     #comgen,
#     ext_prim,
#     ext_sec = repsimintsingle(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);
#
#     SPRICH[i,j] = sprich;
#     RICH[i,j] = rich;
#
#   end #end ltn
# end #end ln
#
# save("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prn_nt.jld","SPRICH",SPRICH,"RICH",RICH);
#
#
# #Load library
# d = load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prn_nt.jld");
# SPRICH = d["SPRICH"];
# RICH = d["RICH"];
#
# richss = Array{Float64}(ln,ltn);
# richsd = Array{Float64}(ln,ltn);
# richcv = Array{Float64}(ln,ltn);
# rrss = Array{Float64}(ln,ltn,reps);
# rrsd = Array{Float64}(ln,ltn,reps);
# for i=1:ln
#   for j=1:ltn
#     reprichss = Array{Float64}(reps);
#     reprichsd = Array{Float64}(reps);
#     reprichcv = Array{Float64}(reps);
#     for r=1:reps
#       traj = copy(SPRICH[i,j][tmax-500:tmax,r])
#       reprichss[r] = mean(traj);
#       reprichsd[r] = std(traj);
#       reprichcv[r] = std(traj)/mean(traj);
#
#       rrss[i,j,r] = mean(traj);
#       rrsd[i,j,r] = var(traj);
#
#     end
#     richss[i,j] = mean(reprichss);
#     richsd[i,j] = mean(reprichsd);
#     richcv[i,j] = mean(reprichcv);
#   end
# end
#
#
# #Plotting the steady state mean
# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/"
# R"""
# pdf(paste($namespace,'fig_sen_prn_nt_mean.pdf',sep=''),width=8,height=6)
# library(RColorBrewer)
# cols <- rev(brewer.pal(10,'Spectral'));
# filled.contour(x=($needvec),y=$n_thresh_vec,z=$richss,xlab='Pr(need)',ylab='n_thresh',key.title=title(main = 'mean'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
# dev.off()
# """
#
# #Plotting steady state variance (SD)
# R"""
# pdf(paste($namespace,'fig_sen_prn_nt_sd.pdf',sep=''),width=8,height=6)
# library(RColorBrewer)
# cols <- rev(brewer.pal(10,'Spectral'));
# filled.contour(x=($needvec),y=$n_thresh_vec,z=$richsd,xlab='Pr(need)',ylab='n_thresh',key.title=title(main = 'sd'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
# dev.off()
# """
#
# #Plotting steady state CV
# R"""
# pdf(paste($namespace,'fig_sen_prn_nt_cv.pdf',sep=''),width=8,height=6)
# library(RColorBrewer)
# cols <- rev(brewer.pal(10,'Spectral'));
# filled.contour(x=($needvec),y=$n_thresh_vec,z=$richcv,xlab='Pr(need)',ylab='n_thresh',key.title=title(main = 'cv'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
# dev.off()
# """
#
# R"""
# pdf(paste($namespace,'fig_taylors.pdf',sep=''),width=6,height=5);
# plot($(vec(rrss[:,:,:])),$(vec(rrsd[:,:,:])),xlab='Richness mean',ylab='Richness variance',pch=16,cex=0.5,log='y')
# dev.off()
# """
#
#
#
# quit()
#
# ##############################
# # pr(make) vs. need_threshold
# ##############################
# #
#
# # using Distributions
# # #using Gadfly
# # using RCall
# # using HDF5
# # using JLD
# #
# # @everywhere using Distributions
# # #@everywhere using Gadfly
# # @everywhere using RCall
# # @everywhere using HDF5
# # @everywhere using JLD
# #
# #
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsim.jl")
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsimintsingle.jl")
# #
# #
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
# #
# #
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonizesingle_func.jl")
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
# # @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")
# #
# #
#
#
# #Read-only variables
# #Establish colonization and extinction rates
# rate_col = 1;
# #Establish thresholds
# a_thresh = 0;
#
# tmax = 2000;
# reps=25;
#
# makevec = collect(0.0001:0.0001:0.001);
# n_thresh_vec = collect(0.0:0.05:0.5);
#
# ln = length(makevec);
# ltn = length(n_thresh_vec);
#
# SPRICH = Array(Array{Int64},ln,ltn);
# RICH =  Array(Array{Int64},ln,ltn);
#
# S = 400;
# ppweight = 1/4;
# sim=true;
# par=true;
#
# for i=1:ln
#   for j=1:ltn
#
#     println("i=",i,"/",ln," j=",j,"/",ltn)
#
#     n_thresh = n_thresh_vec[j]
#     trophicload = 2;
#
#     #Establish community template
#     probs = [
#     p_n=0.008, #This is chosen because there is rich behavior here over n_thresh
#     p_a=0.01,
#     p_m= makevec[i], #needvec[i]/num_play, #1/num_play,
#     p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
#     ]
#
#     sprich,
#     rich,
#     conn,
#     #comgen,
#     ext_prim,
#     ext_sec = repsimintsingle(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);
#
#     SPRICH[i,j] = sprich;
#     RICH[i,j] = rich;
#
#
#   end #end ltn
# end #end ln
#
# save("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prm_nt.jld","SPRICH",SPRICH,"RICH",RICH);
#
#
# # #Load library
# # d = load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prm_nt.jld");
# # SPRICH = d["SPRICH"];
# # RICH = d["RICH"];
#
# richss = Array{Float64}(ln,ltn);
# richsd = Array{Float64}(ln,ltn);
# richcv = Array{Float64}(ln,ltn);
# for i=1:ln
#   for j=1:ltn
#     reprichss = Array{Float64}(reps);
#     reprichsd = Array{Float64}(reps);
#     reprichcv = Array{Float64}(reps);
#     for r=1:reps
#       traj = copy(SPRICH[i,j][tmax-500:tmax,r])
#       reprichss[r] = mean(traj);
#       reprichsd[r] = std(traj);
#       reprichcv[r] = std(traj)/mean(traj);
#     end
#     richss[i,j] = mean(reprichss);
#     richsd[i,j] = mean(reprichsd);
#     richcv[i,j] = mean(reprichcv);
#   end
# end
#
#
# #Plotting the steady state mean
# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/"
# R"""
# pdf(paste($namespace,'fig_sen_prm_nt_mean.pdf',sep=''),width=8,height=6)
# library(RColorBrewer)
# cols <- rev(brewer.pal(10,'Spectral'));
# filled.contour(x=($makevec),y=$n_thresh_vec,z=$richss,xlab='Pr(make)',ylab='n_thresh',key.title=title(main = 'mean'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
# dev.off()
# """
#
# #Plotting steady state variance (SD)
# R"""
# pdf(paste($namespace,'fig_sen_prm_nt_sd.pdf',sep=''),width=8,height=6)
# library(RColorBrewer)
# cols <- rev(brewer.pal(10,'Spectral'));
# filled.contour(x=($makevec),y=$n_thresh_vec,z=$richsd,xlab='Pr(make)',ylab='n_thresh',key.title=title(main = 'sd'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
# dev.off()
# """
#
# #Plotting steady state CV
# R"""
# pdf(paste($namespace,'fig_sen_prm_nt_cv.pdf',sep=''),width=8,height=6)
# library(RColorBrewer)
# cols <- rev(brewer.pal(10,'Spectral'));
# filled.contour(x=($makevec),y=$n_thresh_vec,z=$richcv,xlab='Pr(make)',ylab='n_thresh',key.title=title(main = 'cv'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
# dev.off()
# """
#
#
#



#
# #############################
# pr(make) vs. pr(need)
# #############################


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
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsimintsingle.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonizesingle_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")




#Read-only variables
#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0;

tmax = 2000;
reps=25;

makevec = collect(0.0001:0.0001:0.001);
needvec = collect(0.001:0.001:0.01);

ln = length(needvec);
ltn = length(makevec);

SPRICH = Array(Array{Int64},ln,ltn);
RICH =  Array(Array{Int64},ln,ltn);

S = 400;
ppweight = 1/4;
sim=true;
par=true;

for i=1:ln
  for j=1:ltn

    println("i=",i,"/",ln," j=",j,"/",ltn)

    n_thresh = 0.2;
    trophicload = 2;

    #Establish community template
    probs = [
    p_n=needvec[i], #This is chosen because there is rich behavior here over n_thresh
    p_a=0.01,
    p_m= makevec[j], #needvec[i]/num_play, #1/num_play,
    p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
    ]

    sprich,
    rich,
    conn,
    #comgen,
    ext_prim,
    ext_sec = repsimintsingle(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

    SPRICH[i,j] = sprich;
    RICH[i,j] = rich;


  end #end ltn
end #end ln

save("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prn_prm.jld","SPRICH",SPRICH,"RICH",RICH);


# #Load library
# d = load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_divstats/rich_prn_prm.jld");
# SPRICH = d["SPRICH"];
# RICH = d["RICH"];

richss = Array{Float64}(ln,ltn);
richsd = Array{Float64}(ln,ltn);
richcv = Array{Float64}(ln,ltn);
for i=1:ln
  for j=1:ltn
    reprichss = Array{Float64}(reps);
    reprichsd = Array{Float64}(reps);
    reprichcv = Array{Float64}(reps);
    for r=1:reps
      traj = copy(SPRICH[i,j][tmax-500:tmax,r])
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
pdf(paste($namespace,'fig_sen_prn_prm_mean.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral'));
filled.contour(x=($needvec),y=($makevec),z=$richss,xlab='Pr(need)',ylab='Pr(make)',key.title=title(main = 'mean'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
dev.off()
"""

#Plotting steady state variance (SD)
R"""
pdf(paste($namespace,'fig_sen_prn_prm_sd.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral'));
filled.contour(x=($needvec),y=($makevec),z=$richsd,xlab='Pr(need)',ylab='Pr(make)',key.title=title(main = 'sd'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
dev.off()
"""

#Plotting steady state CV
R"""
pdf(paste($namespace,'fig_sen_prn_prm_cv.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols <- rev(brewer.pal(10,'Spectral'));
filled.contour(x=($needvec),y=($makevec),z=$richcv,xlab='Pr(need)',ylab='Pr(make)',key.title=title(main = 'cv'),color.palette = colorRampPalette(cols, space = "Lab",bias=1),nlevels=50)
dev.off()
"""
