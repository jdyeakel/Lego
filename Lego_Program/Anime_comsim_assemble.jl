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
tmax = 405; #Because we are just packing the community, the max time step should be S+error (from build_template_species.jl)

a_thresh = 0;
trophicload = 2;
ppweight = 1/4;
rate_col = 1;

varvec = [0.0001 0.001 0.002];
numsp = Array(Array{Int64},length(varvec));
potcol = Array(Array{Int64},length(varvec));
trophic = Array(Array{Float64},length(varvec));
maxtrophic = Array(Array{Float64},length(varvec));
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
  mtroph,
  maxtroph = reppak(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  potcol[i] = pot_col;
  numsp[i] = num_sp;
  trophic[i] = mtroph;
  maxtrophic[i] = maxtroph;

end

namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/assemble_pm.jld"
save(namespace,"potcol",potcol,"numsp",numsp,"trophic",trophic,"maxtrophic",maxtrophic);


#Across Need Threshold


S = 400;
reps = 100;
tmax = 405;

a_thresh = 0;
trophicload = 2;
ppweight = 1/4;
rate_col = 1;

varvec = [0.2 0.3 0.4];
numsp = Array(Array{Int64},length(varvec));
potcol = Array(Array{Int64},length(varvec));
trophic = Array(Array{Float64},length(varvec));
maxtrophic = Array(Array{Float64},length(varvec));

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
  mtroph,
  maxtroph = reppak(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  potcol[i] = pot_col;
  numsp[i] = num_sp;
  trophic[i] = mtroph;
  maxtrophic[i] = maxtroph;

end
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/assemble_nt.jld"
save(namespace,"potcol",potcol,"numsp",numsp,"trophic",trophic,"maxtrophic",maxtrophic);

# quit()
# 
# R"plot(seq(1,$tmax)/$num_sp[r],$(pot_col[:,r]/num_sp[r]))"
# R"plot(seq(1,$tmax)/$num_sp[r],$(mtroph[:,r]))"
# 
# 
# R"plot($(pot_col[:,r]/num_sp[r]),$(mtroph[:,r]))"

###Plotting###



varvec_pm = [0.0001 0.001 0.002];
varvec_nt = [0.2 0.3 0.4];
tmax = 405;

dpm = load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/assemble_pm.jld");
potcol_pm = dpm["potcol"];
numsp_pm = dpm["numsp"];
trophic_pm = dpm["trophic"];
dnt = load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/assemble_nt.jld");
potcol_nt = dnt["potcol"];
numsp_nt = dnt["numsp"];
trophic_nt = dnt["trophic"];


namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol_troph.pdf"
R"""
library(RColorBrewer)
cols = brewer.pal(length($varvec_pm),'Set1')
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
pdf($namespace,width=12,height=10);
par(mfrow=c(2,2));
pcol = apply($(potcol_pm[1]),1,mean);
pcolsd = apply($(potcol_pm[1]),1,sd);
ns = mean($(numsp_pm[1]));
plot(seq(1,$tmax)/ns,pcol/ns,type='l',xlim=c(0,1),ylim=c(0,0.25),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1],lwd=3)
polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((pcol-pcolsd)/ns,rev((pcol+pcolsd)/ns)),col=cols[1],border=NA)
legend(x=0.8,y=0.25,legend=$varvec_pm,col=cols,pch=16,title='pr(m)' ,cex=1,bty='n');
"""
for i=2:length(varvec_pm)
  R"""
  pcol = apply($(potcol_pm[i]),1,mean);
  pcolsd = apply($(potcol_pm[i]),1,sd);
  ns = mean($(numsp_pm[i]));
  lines(seq(1,$tmax)/ns,pcol/ns,type='l',xlim=c(0,1),ylim=c(0,0.3),col=cols[$i],lwd=3)
  polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((pcol-pcolsd)/ns,rev((pcol+pcolsd)/ns)),col=cols[$i],border=NA)
  """
end


# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol_nt.pdf"
R"""
library(RColorBrewer)
cols = brewer.pal(length($varvec_nt),'Set1')
#cols[1] = cols[3];
#cols[2] = cols[4];
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
#pdf($namespace,width=12,height=6);
#par(mfrow=c(1,2));
pcol = apply($(potcol_nt[1]),1,mean);
pcolsd = apply($(potcol_nt[1]),1,sd);
ns = mean($(numsp_nt[1]));
plot(seq(1,$tmax)/ns,pcol/ns,type='l',xlim=c(0,1),ylim=c(0,0.25),xlab='Proportion filled',ylab='Proportion potential colonizers',col=cols[1],lwd=3)
polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((pcol-pcolsd)/ns,rev((pcol+pcolsd)/ns)),col=cols[1],border=NA)
legend(x=0.8,y=0.25,legend=$varvec_nt,col=cols,pch=16,title='n_t' ,cex=1,bty='n');
"""
for i=2:length(varvec_nt)
  R"""
  pcol = apply($(potcol_nt[i]),1,mean);
  pcolsd = apply($(potcol_nt[i]),1,sd);
  ns = mean($(numsp_nt[i]));
  lines(seq(1,$tmax)/ns,pcol/ns,type='l',xlim=c(0,1),ylim=c(0,0.3),col=cols[$i],lwd=3)
  polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((pcol-pcolsd)/ns,rev((pcol+pcolsd)/ns)),col=cols[$i],border=NA)
  """
end
#R"dev.off()"

####TROPHIC PLOTS


R"""
library(RColorBrewer)
cols = brewer.pal(length($varvec_pm),'Set1')
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
#pdf($namespace,width=12,height=6);
#par(mfrow=c(1,2));
trm = apply($(trophic_pm[1]),1,mean);
trsd = apply($(trophic_pm[1]),1,sd);
ns = mean($(numsp_pm[1]));
plot(seq(1,$tmax)/ns,trm,type='l',xlim=c(0,1),ylim=c(0,10),xlab='Proportion filled',ylab='Trophic level',col=cols[1],lwd=3)
polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((trm-trsd),rev(trm+trsd)),col=cols[1],border=NA)
legend(x=0.8,y=10,legend=$varvec_pm,col=cols,pch=16,title='pr(m)' ,cex=1,bty='n');
"""
for i=2:length(varvec_pm)
  R"""
  trm = apply($(trophic_pm[i]),1,mean);
  trsd = apply($(trophic_pm[i]),1,sd);
  ns = mean($(numsp_pm[i]));
  lines(seq(1,$tmax)/ns,trm,type='l',xlim=c(0,1),ylim=c(0,0.3),col=cols[$i],lwd=3)
  polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((trm-trsd),rev((trm+trsd))),col=cols[$i],border=NA)
  """
end


# namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_potcol_nt.pdf"
R"""
library(RColorBrewer)
cols = brewer.pal(length($varvec_nt),'Set1')
#cols[1] = cols[3];
#cols[2] = cols[4];
for (i in 1:length(cols)) {cols[i] = paste(cols[i],'70',sep='')}
#pdf($namespace,width=12,height=6);
#par(mfrow=c(1,2));
trm = apply($(trophic_nt[1]),1,mean);
trsd = apply($(trophic_nt[1]),1,sd);
ns = mean($(numsp_nt[1]));
plot(seq(1,$tmax)/ns,trm,type='l',xlim=c(0,1),ylim=c(0,10),xlab='Proportion filled',ylab='Trophic level',col=cols[1],lwd=3)
polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((trm-trsd),rev(trm+trsd)),col=cols[1],border=NA)
legend(x=0.8,y=10,legend=$varvec_nt,col=cols,pch=16,title='n_t' ,cex=1,bty='n');
"""
for i=2:length(varvec_nt)
  R"""
  trm = apply($(trophic_nt[i]),1,mean);
  trsd = apply($(trophic_nt[i]),1,sd);
  ns = mean($(numsp_nt[i]));
  lines(seq(1,$tmax)/ns,trm,type='l',xlim=c(0,1),ylim=c(0,0.3),col=cols[$i],lwd=3)
  polygon(x=c(seq(1,$tmax)/ns,rev(seq(1,$tmax)/ns)),y=c((trm-trsd),rev((trm+trsd))),col=cols[$i],border=NA)
  """
end
R"dev.off()"
