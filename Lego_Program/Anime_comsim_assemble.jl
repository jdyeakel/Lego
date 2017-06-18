using Distributions
# using Gadfly
using RCall
using HDF5
using JLD

@everywhere using Distributions
#@everywhere using Gadfly
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/reppak.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/trophicalc2.jl")



S = 400;
reps = 100;
tmax = 1000;

a_thresh = 0;
trophicload = 2;
ppweight = 1/4;
rate_col = 1;

varvec = [0.0001 0.001 0.002];
numsp = Array(Array{Int64},length(varvec));
potcol = Array(Array{Int64},length(varvec));
trophic = Array(Array{Float64},length(varvec));
for i=1:length(varvec)
  println("i=",i)

  n_thresh = 0.2;

  #Establish community template
  probs = [
  p_n=0.004,
  p_a=0.01,
  p_m= varvec[i], #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]

  #int_m,
  sprich,
  rich,
  conn,
  #comgen,
  ext_prim,
  ext_sec,
  pot_col,
  num_sp,
  mtroph = reppak(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  potcol[i] = pot_col;
  numsp[i] = num_sp;
  trophic[i] = mtroph;

end

namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/assemble_pm.jld"
save(namespace,"potcol",potcol,"numsp",numsp,"trophic",trophic);

# quit()


namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol_comb.pdf"
R"""
library(RColorBrewer)
cols = brewer.pal(length($varvec),'Set1')
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
pdf($namespace,width=12,height=6);
par(mfrow=c(1,2));
pcol = apply($(potcol[1]),1,mean);
pcolsd = apply($(potcol[1]),1,sd);
plot(seq(1,$tmax)/$(numsp[1][1]),pcol/$(numsp[1][1]),type='l',xlim=c(0,1),ylim=c(0,0.25),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1],lwd=3)
polygon(x=c(seq(1,$tmax)/$(numsp[1][1]),rev(seq(1,$tmax)/$(numsp[1][1]))),y=c((pcol-pcolsd)/$(numsp[1][1]),rev((pcol+pcolsd)/$(numsp[1][1]))),col=cols[1],border=NA)
legend(x=0.8,y=0.25,legend=$varvec,col=cols,pch=16,title='pr(m)' ,cex=1,bty='n');
"""
for i=2:length(varvec)
  R"""
  pcol = apply($(potcol[i]),1,mean);
  pcolsd = apply($(potcol[i]),1,sd);
  ns = mean($(numsp[i]));
  lines(seq(1,$tmax)/$(numsp[1][1]),pcol/$(numsp[1][1]),type='l',xlim=c(0,1),ylim=c(0,0.3),col=cols[$i],lwd=3)
  polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((pcol-pcolsd)/ns,rev((pcol+pcolsd)/ns)),col=cols[$i],border=NA)
  """
end
# diff1 = Array{Int64}(0);
# diff2 = Array{Int64}(0);
# for r=1:reps
#   append!(diff1,diff(potcol[1][(find(x->x!=0,potcol[1][:,r])),r]));
#   append!(diff2,diff(potcol[2][(find(x->x!=0,potcol[2][:,r])),r]));
# end
#
# R"""
# hist(c($diff1),xlim=c(0.5,max(c($diff1,$diff2))),ylim=c(0,0.4),col=cols[1],xlab='Niche space expansion',main='',freq=F)
# """
# for i=2:length(varvec)
#   diff2 = Array{Int64}(0);
#   for r=1:reps
#     append!(diff2,diff(potcol[i][(find(x->x!=0,potcol[i][:,r])),r]));
#   end
#   R"""
#   hist(c($diff2),xlim=c(0.5,max(c($diff1,$diff2))),col=cols[$i],xlab='Niche space expansion',main='',add=T,freq=F)
#   """
# end
# R"dev.off()"
#
#
#
# quit()

#
# workspace();
# using Distributions
# # using Gadfly
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
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/reppak.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")
#
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
#

#Across Need Threshold


S = 400;
reps = 100;
tmax = 1000;

a_thresh = 0;
trophicload = 2;
ppweight = 1/4;
rate_col = 1;

varvec = [0.2 0.3 0.4];
numsp = Array(Array{Int64},length(varvec));
potcol = Array(Array{Int64},length(varvec));
trophic = Array(Array{Float64},length(varvec));
for i=1:length(varvec)
  println("i=",i)

  n_thresh = varvec[i];

  #Establish community template
  probs = [
  p_n=0.004,
  p_a=0.01,
  p_m= 0.001, #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]

  #int_m,
  sprich,
  rich,
  conn,
  #comgen,
  ext_prim,
  ext_sec,
  pot_col,
  num_sp,
  mtroph = reppak(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  potcol[i] = pot_col;
  numsp[i] = num_sp;
  trophic[i] = mtroph;

end
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/assemble_nt.jld"
save(namespace,"potcol",potcol,"numsp",numsp,"trophic",trophic);

# quit()


# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol_nt.pdf"
R"""
library(RColorBrewer)
cols = brewer.pal(length($varvec),'Set1')
#cols[1] = cols[3];
#cols[2] = cols[4];
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
#pdf($namespace,width=12,height=6);
#par(mfrow=c(1,2));
pcol = apply($(potcol[1]),1,mean);
pcolsd = apply($(potcol[1]),1,sd);
plot(seq(1,$tmax)/$(numsp[1][1]),pcol/$(numsp[1][1]),type='l',xlim=c(0,1),ylim=c(0,0.25),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1],lwd=3)
polygon(x=c(seq(1,$tmax)/$(numsp[1][1]),rev(seq(1,$tmax)/$(numsp[1][1]))),y=c((pcol-pcolsd)/$(numsp[1][1]),rev((pcol+pcolsd)/$(numsp[1][1]))),col=cols[1],border=NA)
legend(x=0.8,y=0.25,legend=$varvec,col=cols,pch=16,title='n_t' ,cex=1,bty='n');
"""
for i=2:length(varvec)
  R"""
  pcol = apply($(potcol[i]),1,mean);
  pcolsd = apply($(potcol[i]),1,sd);
  ns = mean($(numsp[i]));
  lines(seq(1,$tmax)/$(numsp[1][1]),pcol/$(numsp[1][1]),type='l',xlim=c(0,1),ylim=c(0,0.3),col=cols[$i],lwd=3)
  polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((pcol-pcolsd)/ns,rev((pcol+pcolsd)/ns)),col=cols[$i],border=NA)
  """
end
# diff1 = Array{Int64}(0);
# diff2 = Array{Int64}(0);
# for r=1:reps
#   append!(diff1,diff(potcol[1][(find(x->x!=0,potcol[1][:,r])),r]));
#   append!(diff2,diff(potcol[2][(find(x->x!=0,potcol[2][:,r])),r]));
# end
#
# R"""
# hist(c($diff1),xlim=c(0.5,max(c($diff1,$diff2))),ylim=c(0,0.4),col=cols[1],xlab='Niche space expansion',main='',freq=F)
# """
# for i=2:length(varvec)
#   diff2 = Array{Int64}(0);
#   for r=1:reps
#     append!(diff2,diff(potcol[i][(find(x->x!=0,potcol[i][:,r])),r]));
#   end
#   R"""
#   hist(c($diff2),xlim=c(0.5,max(c($diff1,$diff2))),col=cols[$i],xlab='Niche space expansion',main='',add=T,freq=F)
#   """
# end
R"dev.off()"



#
#
# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol_nt.pdf"
# R"""
# library(RColorBrewer)
# cols = brewer.pal(max(5,length($varvec)),'Set1')
# #cols[1] = cols[3];
# #cols[2] = cols[4];
# for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
# #pdf($namespace,width=12,height=6);
# par(mfrow=c(1,2));
# pcol = $(potcol[1]);
# plot(seq(1,$tmax)/$(numsp[1][1]),pcol[,1]/$(numsp[1][1]),type='l',xlim=c(0,1),ylim=c(0,0.3),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1])
# for (r in 2:$reps) {
#   ns = $(numsp[1]);
#   lines(seq(1,$tmax)/ns[r],pcol[,r]/ns[r],col=cols[1])
# }
# legend(x=0.8,y=0.3,legend=$varvec,col=cols,pch=16,title='need thresh',cex=1,bty='n');
# """
# for i=2:length(varvec)
#   R"""
#   pcol = $(potcol[i]);
#   ns = $(numsp[i]);
#   for (r in 1:$reps) {
#     lines(seq(1,$tmax)/ns[r],pcol[,r]/ns[r],col=cols[$i])
#   }
#   """
# end
# diff1 = Array{Int64}(0);
# diff2 = Array{Int64}(0);
# for r=1:reps
#   append!(diff1,diff(potcol[1][(find(x->x!=0,potcol[1][:,r])),r]));
#   append!(diff2,diff(potcol[2][(find(x->x!=0,potcol[2][:,r])),r]));
# end
# R"""
# hist(c($diff1),xlim=c(0.5,max(c($diff1,$diff2))),ylim=c(0,0.2),col=cols[1],xlab='Niche space expansion',breaks=10,main='',freq=F)
# hist(c($diff2),xlim=c(0.5,max(c($diff1,$diff2))),col=cols[2],xlab='Niche space expansion',breaks=10,main='',add=T,freq=F)
# #dev.off()
# """
#
# #polygon(c(seq(1,$tmax)/$(numsp[1][1]),rev(seq(1,$tmax)/$(numsp[1][1]))),c((pcol-pcolsd)/$(numsp[1][1]),rev((pcol+pcolsd)/$(numsp[1][1]))),col=col[1])
# #lines(seq(1,$tmax)/$(numsp[1][1]),(pcol+pcolsd)/$(numsp[1][1]),col=cols[1])
# #lines(seq(1,$tmax)/$(numsp[1][1]),(pcol-pcolsd)/$(numsp[1][1]),col=cols[1])
