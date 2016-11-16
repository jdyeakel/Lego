
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
n_thresh = 0.2;
tmax = 800;
reps=500;

needvec = collect(0.5:1.0:5.5);
ln = length(needvec);

SPARSECG = Array(Array{Float64},ln);
RvsE = Array(Array{Float64},ln);
RvsEsd = Array(Array{Float64},ln);
rich_ss = zeros(Float64,ln);
COMGEN = Array(Array{Int64},ln);
INT =  Array(Array{Char},ln);
num_play = 500;
ppweight = 1/4;
sim=true;
par=true;
#Loop over need value
for i=1:ln


  #Establish community template
  probs = [
  p_n=0.003,
  p_a=0.01,
  p_m=needvec[i]/num_play,  #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]

  int_m, sprich, rich, conn, comgen, ext_prim, ext_sec = repsim(num_play,reps,tmax,a_thresh,n_thresh,rate_col,probs,ppweight,sim,par);

  COMGEN[i] = comgen;
  INT[i] = int_m;

  #Number of extinctions for a community of a given species richness
  mext = zeros(maximum(sprich),2);
  sdext = zeros(maximum(sprich),2);
  for j=1:maximum(sprich)
    numsp=find(x->x==j,sprich); #find the location of sprichness=j on the matrix
    lnumsp = length(numsp);
    randsamp = sample(collect(1:lnumsp),min(lnumsp,300),replace=true);
    mext[j,1] = j;
    mext[j,2] = mean(ext_sec[numsp][randsamp]+ext_prim[numsp][randsamp]);
    sdext[j,1] = j;
    sdext[j,2] = std(ext_sec[numsp][randsamp]+ext_prim[numsp][randsamp]);
  end
  RvsE[i] = mext;
  RvsEsd[i] = sdext;
  rich_ss[i] = mean(sprich[(tmax-100):(tmax),:]);
end





################
# SAVING
################

#SAVE DATA FOR NEED SIMULATIONS
#Compression and storage of COMGEN
for i=1:ln
  for r=1:reps
    sparsecg = sparse(COMGEN[i][r,:,:]);
    namespace = string("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/comgen_n/sparsecg_i",i,"_r",r,".jld");
    save(namespace,"sm",sparsecg);
  end
end
#The rest of the data
save("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/overn.jld","INT",INT,"rve",RvsE,"rvesd",RvsEsd,"rich",rich_ss);


#SAVE DATA FOR MAKE SIMULATIONS
#Compression and storage of COMGEN
for i=1:ln
  for r=1:reps
    sparsecg = sparse(COMGEN[i][r,:,:]);
    namespace = string("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/comgen_m/sparsecg_i",i,"_r",r,".jld");
    save(namespace,"sm",sparsecg);
  end
end
#Save Data
save("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/overm.jld","INT",INT,"rve",RvsE,"rvesd",RvsEsd,"rich",rich_ss);




################
# LOADING
################

#LOAD NEED DATA
COMGEN = Array(Array{Int64},ln);
for i=1:ln
  comgen = Array{Int64}(reps,num_play,tmax);
  for r=1:reps
    namespace = string("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/comgen_n/sparsecg_i",i,"_r",r,".jld");
    libr=load(namespace);
    comgen[r,:,:] = libr["sm"];
  end
  COMGEN[i] = comgen;
end
d = load("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/overn.jld");
#This loads the dictionary
RvsE = d["rve"];
rich_ss = d["rich"];


#LOAD MAKE DATA
COMGEN = Array(Array{Int64},ln);
for i=1:ln
  comgen = Array{Int64}(reps,num_play,tmax);
  for r=1:reps
    namespace = string("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/comgen_m/sparsecg_i",i,"_r",r,".jld");
    libr=load(namespace);
    comgen[r,:,:] = libr["sm"];
  end
  COMGEN[i] = comgen;
end
d = load("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_ext/overm.jld");
#This loads the dictionary
RvsE = d["rve"];
rich_ss = d["rich"];








################
# PLOTTING NEED DATA
################


#Plot Extinctions vs. Richness
R"""
library(RColorBrewer);
cols <- brewer.pal(10,'Spectral');
x=$(RvsE[1][:,1]);
y=$(RvsE[1][:,2]);
pdf('/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_rve.pdf',width=5,height=4)
plot(x,y,type='l',ylim=c(0,10),xlim=c(0,70),col=cols[1],lwd=2,
  xlab="Species richness",ylab="Mean number of extinctions");
points($(rich_ss[1]),0,pch=16,col=cols[1]);
segments(x0=$(rich_ss[1]),y0=0,x1=$(rich_ss[1]),y1=-1,col=cols[1]);
legvec = seq(0.5,5,0.5)/500;
legend(x=60,y=10,legend=legvec,col=cols,pch=16,title='Prob(need)',cex=0.5,bty='n');
"""
for i=2:10
  x=(RvsE[i][:,1]);
  y=(RvsE[i][:,2]);
  richss=rich_ss[i];
  R"lines($x,$y,col=cols[$i],lwd=2);
  points($richss,0,pch=16,col=cols[$i]);
  segments(x0=$richss,y0=0,x1=$richss,y1=-1,col=cols[$i]);
  "
end
R"dev.off()"

#Plot Extinctions vs. Richness Normalized to 1
R"""
library(RColorBrewer);
cols <- brewer.pal(10,'Spectral');
x=$(RvsE[1][:,1]);
y=$(RvsE[1][:,2])/max($(RvsE[1][:,2]));
pdf('/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_rvenorm.pdf',width=5,height=4)
plot(x,y,type='l',ylim=c(0,1),xlim=c(0,80),col=cols[1],lwd=2,
  xlab="Species richness",ylab="Normalized mean number of extinctions");
points($(rich_ss[1]),0,pch=16,col=cols[1]);
segments(x0=$(rich_ss[1]),y0=0,x1=$(rich_ss[1]),y1=-1,col=cols[1]);
legvec = seq(0.5,5,0.5)/500;
legend(x=70,y=1,legend=legvec,col=cols,pch=16,title='Prob(need)',cex=0.5,bty='n');
"""
for i=2:10
  x=(RvsE[i][:,1]);
  y=(RvsE[i][:,2])/maximum((RvsE[i][:,2]));
  richss=rich_ss[i];
  R"lines($x,$y,col=cols[$i],lwd=2);
  points($richss,0,pch=16,col=cols[$i]);
  segments(x0=$richss,y0=0,x1=$richss,y1=-1,col=cols[$i]);
  "
end
R"dev.off()"



################
# PLOTTING MAKE DATA
################


#Plot Extinctions vs. Richness
R"""
library(RColorBrewer);
cols <- brewer.pal(10,'Spectral');
x=$(RvsE[1][:,1]);
y=$(RvsE[1][:,2]);
#pdf('/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_rveM.pdf',width=5,height=4)
plot(x,y,type='l',ylim=c(0,10),xlim=c(0,100),col=cols[1],lwd=2,
  xlab="Species richness",ylab="Mean number of extinctions");
points($(rich_ss[1]),0,pch=16,col=cols[1]);
segments(x0=$(rich_ss[1]),y0=0,x1=$(rich_ss[1]),y1=-1,col=cols[1]);
legvec = seq(0.5,5,0.5)/500;
legend(x=90,y=10,legend=legvec,col=cols,pch=16,title='Prob(make)',cex=0.5,bty='n');
"""
for i=2:10
  x=(RvsE[i][:,1]);
  y=(RvsE[i][:,2]);
  richss=rich_ss[i];
  R"lines($x,$y,col=cols[$i],lwd=2);
  points($richss,0,pch=16,col=cols[$i]);
  segments(x0=$richss,y0=0,x1=$richss,y1=-1,col=cols[$i]);
  "
end
R"dev.off()"

#Plot Extinctions vs. Richness Normalized to 1
R"""
library(RColorBrewer);
cols <- brewer.pal(10,'Spectral');
x=$(RvsE[1][:,1]);
y=$(RvsE[1][:,2])/max($(RvsE[1][:,2]));
pdf('/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_rvenormM.pdf',width=5,height=4)
plot(x,y,type='l',ylim=c(0,1),xlim=c(0,110),col=cols[1],lwd=2,
  xlab="Species richness",ylab="Normalized mean number of extinctions");
points($(rich_ss[1]),0,pch=16,col=cols[1]);
segments(x0=$(rich_ss[1]),y0=0,x1=$(rich_ss[1]),y1=-1,col=cols[1]);
legvec = seq(0.5,5,0.5)/500;
legend(x=95,y=1,legend=legvec,col=cols,pch=16,title='Prob(make)',cex=0.5,bty='n');
"""
for i=2:10
  x=(RvsE[i][:,1]);
  y=(RvsE[i][:,2])/maximum((RvsE[i][:,2]));
  richss=rich_ss[i];
  R"lines($x,$y,col=cols[$i],lwd=2);
  points($richss,0,pch=16,col=cols[$i]);
  segments(x0=$richss,y0=0,x1=$richss,y1=-1,col=cols[$i]);
  "
end
R"dev.off()"





richextplot = plot(
[layer(x=RvsE[j][:,1],y=RvsE[j][:,2], Geom.line, Theme(default_color=colorant"black",default_point_size=2pt, highlight_width = 0pt)) for j in 1:ln]...,
Guide.xlabel("Species richness"),Guide.ylabel("Mean total extinctions"));
