using Distributions
using Gadfly
using RCall
#using PyCall
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")

#Establish community template
num_play = 100;
probs = [
p_n=0.1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
sim=true;
ppweight=1/3;
int_m, sp_m, t_m, tp_m, tind_m, mp_m, simvalue = build_template_degrees(num_play,probs,ppweight,sim);

#Establish colonization and extinction rates
rate_col = 0.2;
rate_ext = 0.1;

#Establish thresholds
a_thresh = 0.2;
n_thresh = 0.4;

rep = 100;
CID = (Array{Int64,1})[];
rich = Array{Int64}(rep);
for r = 1:rep
  #The add-until-full simulation
  #Creating a new int_m each time
  sim=true;
  ppweight=1/3;
  int_m, sp_m, t_m, tp_m, tind_m, mp_m, simvalue = build_template_degrees(num_play,probs,ppweight,sim);
  cid, c_m, crev_m, com_tp, com_tind, com_mp = initiate_comm_func(int_m,tp_m,tind_m,mp_m);
  status = "open";
  while status == "open"
    status,cid,c_m,crev_m,com_tp,com_tind,com_mp = colonize_func(int_m,tp_m,tind_m,mp_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp);
  end
  # length(unique(cid))-length(cid)
  rich[r] = length(cid);
  println("Richness = ",rich[r])
  push!(CID,copy(cid));
end
csum = (Array{Int64,1})[];
for i=1:rep
  push!(csum,cumsum(CID[i]));
end
#Visualize the assembly process
assembplot = plot(
[layer(y=csum[j],x=collect(1:length(csum[j])), Geom.line, Theme(default_color=colorant"black")) for j in 1:rep]...,
Guide.xlabel("Steps"),Guide.ylabel("Summed species ID"));

draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_fullassemb.pdf", 5inch, 4inch), assembplot)



############################
############################
# COLONIZATION + EXTINCTION
############################
############################
using Distributions
#using Gadfly
using RCall
#using PyCall
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")


#Establish community template
num_play = 500;
probs = [
p_n=0.010,
p_a=0.01,
p_m=1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
#int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,probs);

#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0.0;
n_thresh = 0.1;
trophicload = 2;
tmax = 3000;
CID = (Array{Int64,1})[];
fwt = Array(Array{Int64},tmax);
com_probs = Array{Float64}(tmax,4);
rich = Array{Int64}(tmax);
sprich = Array{Int64}(tmax);
conn = Array{Float64}(tmax);
ext_prim = Array{Int64}(tmax);
ext_sec = Array{Int64}(tmax);
comgen =zeros(Int64,tmax,num_play);
ppweight = 1/4;
sim=false;
par=false;

window = 20;
mrich = zeros(0);
colrate = zeros(0);
extrate = zeros(0);
tictoc=0;
sumcol=0;
sumext=0;

@time int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_degrees(num_play,probs,ppweight,sim,par);
cid, c_m, crev_m, com_tp, com_tind, com_mp, com_mind = initiate_comm_func(int_m,tp_m,tind_m,mp_m,mind_m);
status = "open"; #I don't think we need this now
@time for t = 1:tmax
  #The add-until-full simulation
  #Creating a new int_m each time

  #Colonize with some probability
  rcol = rand();
  if rcol < rate_col && status == "open"
    status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind = colonize_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
    if status == "open"
      sumcol=sumcol+1;
    end
  end

  #Always run extinction code because probabilities are assessed within
  status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,extinctions = extinct_func2(int_m,tp_m,a_thresh,n_thresh,trophicload,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
  #status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,extinctions = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,simvalue);
	#Save primary and secondary extinction information

  ext_prim[t] = extinctions[1];
	ext_sec[t] = extinctions[2];
  sumext = sumext + sum(extinctions);

  S = length(spcid);
  # length(unique(cid))-length(cid)
  conn[t] = (sum(com_tp))/(S^2);


  rich[t] = length(cid);
  sprich[t] = length(spcid);

  #println("Richness = ",rich[t])

  push!(CID,copy(cid));
  comgen[t,cid] = 1;
  fwt[t] = com_tind;

  # int_mc = copy(int_m);
  # for i=1:num_play
  #   int_mc[i,i]='0';
  # end


  if mod(t,window)==0
    tictoc=tictoc+1;
    mt = mean(sprich[(t-window+1):t]);
    append!(mrich,mt);
    append!(colrate,sumcol/mt);
    append!(extrate,sumext/mt);
    #Reset
    sumcol=0;
    sumext=0;
  end



  # #Counting probabilities of a, n, i, m, e within simulated communities
  # com_int = sum(vec(int_mc[cid,cid]).==vec(['a', 'n', 'i', 'm'])',1)./(length(cid)^2);
  # temp_int = sum(vec(int_mc).==vec(['a', 'n', 'i', 'm'])',1)./(length(int_mc));
  # com_probs[t,:] = (com_int./temp_int) #/ sum((com_int./temp_int));
  #
end
#writedlm("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen.csv",comgen);


#namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/"


R"""
maxsp = max($rich);
plot($sprich,type='l',xlab='Time',ylab='Species diversity',ylim=c(0,maxsp))
lines($rich,type='l',lty=2,cex.lab=1.5,cex.axis=1.3)
recol=$(find(x->x==0,sprich));
points(recol,rep(0,length(recol)))
"""



R"""
#pdf(paste($namespace,"fig_colextrates.pdf",sep=""),width=6,height=5)
library(RColorBrewer)
cols <- brewer.pal(3,"Set1");
maxvalue = max(c($colrate,$extrate))
plot($mrich,$colrate,pch=16,type='o',cex=0.5,col=cols[1],xlab='Species richness',ylab='Rate',ylim=c(0,maxvalue))
points($mrich,$extrate,pch=16,type='o',cex=0.5,col=cols[2])
meansprich=mean($(sprich[Int(round(tmax/2)):tmax]));
segments(x0=meansprich,y0=-1,x1=meansprich,y1=maxvalue+1,lty=2)
#dev.off()
"""

# 
# c=sort(CID[tmax]);
# dups = Array{Int64}(0);
# for i=2:length(c)
#   if c[i]-c[i-1]==0
#     push!(dups,c[i])
#   end
# end



R"""
plot($conn,type='l',log='y',ylab='Trophic connectance')
"""

R"""
plot(($rich-$sprich),$conn,xlab='Number of objects',ylab='Connectance')
"""

R"""
plot(($rich-$sprich),$sprich,xlab='Number of objects',ylab='Species richness')
"""

R"""
plot($ext_prim,$ext_sec)
"""

R"""
plot($(com_probs[:,1]),type='l',col='red',ylim=c(0,15))
lines($(com_probs[:,2]),col='blue')
lines($(com_probs[:,3]),col='black')
lines($(com_probs[:,4]),col='purple')
"""

R"""
library(igraph)
library(RColorBrewer)
pal <- brewer.pal(3,"Set1")
fw <- as.matrix($tp_m)
#Eliminate zeros
todel<-which((apply(fw,2,sum)*apply(fw,1,sum))==0)
fw <- fw[-todel,-todel];
fw_g <- graph.adjacency(as.matrix(fw))
basal_pos <- 1
num_play = dim(fw)[1]
trophic <- sapply(1:vcount(fw_g),function(x){mean(shortest.paths(fw_g,basal_pos,which(fw[x,]==1)))+0})
#trophic[which(trophic==Inf)] <- 0
trophic[which(trophic=="NaN")] <- 0
coords <- cbind(runif(vcount(fw_g)),trophic); coords[basal_pos,] <- c(0.5,trophic[basal_pos])
par(mar=c(1,1,1,1))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.5,
     main=ecount(fw_g)/(num_play^2),vertex.label=NA,
     vertex.color=pal[2]) #,layout=coords
"""



csum = (Array{Int64,1})[];
for i=1:rep
  push!(csum,cumsum(CID[i]));
end
#Diversity through time
divplot = plot(y=rich,x=collect(1:tmax), Geom.line, Theme(default_color=colorant"black"),
Guide.xlabel("Time"),Guide.ylabel("Species richness"));
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_diversity.pdf", 5inch, 4inch), divplot)

#Food web connectance through time
conplot = plot(y=conn,x=collect(1:tmax), Geom.line, Scale.x_log10,Scale.y_log10, Theme(default_color=colorant"black"),
Guide.xlabel("Time"),Guide.ylabel("Food web connectance"));
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_connectance.pdf", 5inch, 4inch), divplot)

plot(x=rich,y=conn,Geom.line,Scale.x_log10,Scale.y_log10,Geom.point)




####################################
####################################
# PARALLEL COLONIZATION + EXTINCTION
####################################
####################################

using Distributions
using Gadfly
using RCall

@everywhere using Distributions
@everywhere using Gadfly
@everywhere using RCall

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
tmax = 500;
reps=100;

num_play = 1000;
#Shared variables
sprich = SharedArray{Int64}(tmax,reps);
rich = SharedArray{Int64}(tmax,reps);
conn = SharedArray{Float64}(tmax,reps);
comgen = SharedArray{Int64}(reps,num_play,tmax);
ext_prim = SharedArray{Int64}(tmax,reps);
ext_sec = SharedArray{Int64}(tmax,reps);

#Establish community template
probs = [
p_n=1/num_play,
p_a=0.01,
p_m=1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
ppweight = 1/4;
sim=true;
par=true;

@time int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_degrees(num_play,probs,ppweight,sim,par);

@sync @parallel for r=1:reps
  #Establish community template
  cid, c_m, crev_m, com_tp, com_tind, com_mp, com_mind = initiate_comm_func(int_m,tp_m,tind_m,mp_m,mind_m);
  for t=1:tmax
    status = "open";
    #Colonize with some probability
    rcol = rand();
    if rcol < rate_col && status == "open"
      status,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind = colonize_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
    end
    #Always run extinction module because probabilities are assessed within
    status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,extinctions = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,simvalue);
    #Save primary and secondary extinction information
  	ext_prim[t,r] = extinctions[1];
  	ext_sec[t,r] = extinctions[2];
    sprich[t,r] = length(spcid);
    # length(unique(cid))-length(cid)
    conn[t,r] = (sum(com_tp))/(sprich[t,r]^2);
    rich[t,r] = length(cid);
    comgen[r,cid,t] = 1;
  end
end

#Similarity posthoc analysis
jacc = SharedArray{Float64}(reps,reps,tmax);
meansim = SharedArray{Float64}(reps,tmax);
maxsim = SharedArray{Float64}(reps,tmax);
@sync @parallel for t=1:tmax
  commat = copy(comgen[:,:,t]);
  for i=1:reps
    x = copy(commat[i,:]);
    for j=1:reps
      y = copy(commat[j,:]);
      if j >= i
        # sim = dot(mati,matj)/mag(i)
        jacc[i,j,t] = 1 - sum(min(x, y)) / sum(max(x, y));
        jacc[j,i,t] = copy(jacc[i,j,t]);
      end
    end
    meansim[i,t] = mean(1-jacc[i,:,t]);
    maxsim[i,t] = sort(1-jacc[i,:,t],rev=true)[2];
  end
  print(t)
end


#Save the data
for t=1:tmax
  namespace = string("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen_t",t,".csv");
  writedlm(namespace,comgen[:,:,t]);
end

##########
#ANALYSIS#
##########

#NOTE: rewrite plots in R?

#Visualize richness over time
# tmax = 500
richplot = plot(
[layer(y=sprich[1:tmax,j],x=collect(1:tmax), Geom.line, Theme(default_color=colorant"gray")) for j in 1:reps]...,
Guide.xlabel("Time"),Guide.ylabel("Richness"));
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_diversity.pdf", 5inch, 4inch), richplot)

#Connectance Plot
connplot = plot(
[layer(y=conn[:,j],x=collect(1:tmax), Geom.line, Theme(default_color=colorant"gray")) for j in 1:reps]...,
Guide.xlabel("Time"),Guide.ylabel("Connectance"),Scale.x_log10,Scale.y_log10);
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_connectance.pdf", 5inch, 4inch), connplot)

#Similarity plot
simplot = plot(
[layer(y=maxsim[j,:],x=collect(1:tmax), Geom.line, Theme(default_color=colorant"gray")) for j in 1:reps]...,
Guide.xlabel("Time"),Guide.ylabel("Similarity"),Scale.x_log10)


#Similarity plot
simmaxplot = plot(
[layer(y=meansim[j,:],x=maxsim[j,:], Geom.point, Theme(default_color=colorant"black",default_point_size=1pt, highlight_width = 0pt)) for j in 1:reps]...,
Guide.xlabel("Maximum similarity"),Guide.ylabel("Mean similarity"));
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_simmax.pdf", 5inch, 4inch), simmaxplot)

#Extinction plot with jitter
extplot =  plot(x=vec(ext_prim).+rand(collect(0:0.001:0.4),tmax*reps),y=vec(ext_sec).+rand(collect(0:0.001:0.4),tmax*reps), Geom.point, Theme(default_color=colorant"black",default_point_size=0.5pt, highlight_width = 0pt),
Guide.xlabel("Num primary extinctions"),Guide.ylabel("Num secondary extinctions"))
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_ext.pdf", 5inch, 4inch), extplot)

maxprim = maximum(ext_prim);
msec = zeros(maxprim+1)
sdsec = zeros(maxprim+1)
for i=0:maxprim
  msec[i+1] = mean(ext_sec[find(x->x==i,ext_prim)]);
  sdsec[i+1] = var(ext_sec[find(x->x==i,ext_prim)]);
end
plot(x=collect(0:maxprim),y=sdsec)

plot(x=msec,y=sdsec)

plot(x=vec(sprich).+rand(collect(0:0.001:0.4),tmax*reps),y=(vec(ext_prim).+vec(ext_sec)).+rand(collect(0:0.001:0.4),tmax*reps), Geom.point, Theme(default_color=colorant"black",default_point_size=0.5pt, highlight_width = 0pt))


#Number of extinctions for a community of a given species richness
mext = zeros(maximum(sprich));
for i=1:maximum(sprich)
  numsp=find(x->x==i,sprich);
  lnumsp = length(numsp);
  randsamp = sample(collect(1:lnumsp),100,replace=true)
  mext[i] = mean(ext_sec[numsp][randsamp]+ext_prim[numsp][randsamp]);
end
richextplot = plot(x=collect(1:maximum(sprich)),y=mext, Geom.point, Theme(default_color=colorant"black",default_point_size=2pt, highlight_width = 0pt),
Guide.xlabel("Species richness"),Guide.ylabel("Mean total extinctions"));
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_richext.pdf", 5inch, 4inch), richextplot)
