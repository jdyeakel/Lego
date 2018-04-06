############################
############################
# COLONIZATION + EXTINCTION
############################
############################
loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/loadfuncs.jl");

#Establish community template
S = 400;
# S = 400;
probs = [
# p_n=0.04,
# p_a=0.01,
# p_m=0.04,
p_n=0.004,
p_a=0.01,
p_m=0.002,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
#int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,probs);

#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0.0;
n_thresh = 0.2;
trophicload = 2;
tmax = 5000;


CID = (Array{Int64,1})[];
SPCID = (Array{Int64,1})[];
fwt = Array{Array{Int64}}(tmax);
fwtind = Array{Array{Int64}}(tmax);
fwm = Array{Array{Int64}}(tmax);
fwmind = Array{Array{Int64}}(tmax);
com_probs = Array{Float64}(tmax,4);
rich = Array{Int64}(tmax);
sprich = Array{Int64}(tmax);
conn = Array{Float64}(tmax);
connind = Array{Float64}(tmax);
ext_prim = Array{Int64}(tmax);
ext_sec = Array{Int64}(tmax);
pot_col = Array{Array{Int64}}(tmax);
tl = Array{Float64}(tmax);
tw = Array{Float64}(tmax);
twsd = Array{Float64}(tmax);
twind = Array{Float64}(tmax);
twindsd = Array{Float64}(tmax);


ppweight = 1/4;
sim=false;
par=false;
calcpotcol = false;

window = 20;
mrich = zeros(0);
colrate = zeros(0);
extrate = zeros(0);
tictoc=0;
sumcol=0;
sumext=0;

@time int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_species(S,probs,ppweight);

num_play = length(diag(int_m));
comgen =zeros(Int64,tmax,num_play);
cid = Array{Int64}(0);

cid, c_m, crev_m, com_tp, com_tind, com_mp, com_mind = initiate_comm_func(int_m,tp_m,tind_m,mp_m,mind_m);
status = "open"; #I don't think we need this now
@time for t = 1:tmax
  #Print every 1000 timesteps
  if mod(t,1000)==0
    println(string("t=",t))
  end
  #The add-until-full simulation
  #Creating a new int_m each time

  #Colonize with some probability
  rcol = rand();
  if rcol < rate_col && status == "open"

    #Run the faster colonization function if calcpotcol is set to false

    if calcpotcol == true
      status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,potcol = colonize_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
    else
      status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind= colonizesingle_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
    end

    if status == "open"
      sumcol=sumcol+1;
    end
  end

  #Always run extinction code because probabilities are assessed within
  status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,extinctions = extinct_func2(int_m,tp_m,a_thresh,n_thresh,trophicload,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
  #status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,extinctions = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,simvalue);
	#Save primary and secondary extinction information

  rich[t] = length(cid);
  sprich[t] = length(spcid);


end

R"""
plot($(collect(1:tmax)),$rich,type='l',lty=2)
lines($(collect(1:tmax)),$sprich)
"""


#writedlm("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen.csv",comgen);
obrich = rich.-sprich;
ext_tot = ext_prim .+ ext_sec;

#namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/"
if calcpotcol==true
  lpotcol=zeros(Int64,tmax)
  for i=1:tmax
    lpotcol[i]=length(pot_col[i])
  end
end


##########################
#Proportion of stable webs
##########################
#Turn on parallelization
par=false;
psw = zeros(Float64,tmax);
pswind = zeros(Float64,tmax);
for t=1:tmax
    sigma = 0.4;
    reps = 10;
    psw[t] = PSWebs(fwt[t],fwm[t],sigma,reps,par);
    pswind[t] = PSWebs(fwtind[t],fwmind[t],sigma,reps,par);
    if mod(t,100)==0
      println(string("t=",t))
    end
end

sscomm = mean(sprich[Int64(floor(tmax/2))]);
sscomm_ob = mean(rich[Int64(floor(tmax/2))]);
R"""
library(RColorBrewer)
pal = brewer.pal(3,"Set1")
par(mfrow=c(2,1))
plot($(movingaverage(psw,10)),col=pal[1],type='l',ylim=c(0,1.2))
lines($(movingaverage(pswind,10)),col=pal[2])
lines($(movingaverage(sprich/sscomm,10)),col='black')
lines($(movingaverage(rich/sscomm,10)),col='black',lty=2)
plot($sprich,$psw,col=pal[1],pch=16,cex=0.8)
points($sprich,$pswind,col=pal[2],pch=16,cex=0.8)
"""

#PSW at time t vs. Primary extinctions at t+1
delay=1;
R"""
plot($(pswind[1:tmax-delay]),$(ext_prim[(delay+1):tmax]),ylab='Number prim extinctions @ time t+1',xlab='PSW @ time t',col=pal[1],pch=16,cex=0.8)
#points(jitter($(pswind[1:tmax-delay])),$(ext_prim[delay+1:tmax]),col=pal[2],pch=16,cex=0.8)*/
"""

delay=1;
R"""
plot(jitter($(pswind[1:tmax-delay])),$(ext_prim[delay+1:tmax]),ylab=paste('Number prim extinctions @ time t+',$delay,sep=''),xlab='PSW @ time t')
points(jitter($(pswind[1:tmax-delay])),$(ext_sec[delay+1:tmax]),ylab=paste('Number prim extinctions @ time t+',$delay,sep=''),xlab='PSW @ time t',pch=16,cex=0.8)
"""

delay=1;
R"plot(jitter($(pswind[1:tmax-delay])),$(ext_sec[delay+1:tmax])/$(sprich[delay+1:tmax]),col='blue',ylab=paste('Number prim extinctions @ time t+',$delay,sep=''),xlab='PSW @ time t')"


namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/";
R"""
pdf(paste($namespace,"fig_trophicwidth.pdf",sep=""),width=8,height=4)
par(mfrow=c(1,2));plot($tw,type='l',log='x',xlab='Time',ylab='Trophic overlap mean')
plot($twsd,type='l',log='x',xlab='Time',ylab='Trophic overlap SD')
dev.off()
"""


R"plot($twind,type='l',log='x',xlab='Time',ylab='Indirect trophic width')"
R"plot($twindsd,type='l',log='x',xlab='Time',ylab='Indirect trophic width SD')"



#Calculate the final distribution of trophic levels
tl = trophicalc(com_tp);


namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/";
R"""
#pdf(paste($namespace,"fig_traj.pdf",sep=""),width=6,height=5)*/
maxsp = max($rich);
plot($sprich,type='l',xlab='Time',ylab='Species diversity',ylim=c(0,maxsp),log='x')
lines($rich,type='l',lty=2,cex.lab=1.5,cex.axis=1.3)
recol=$(find(x->x==0,sprich));
points(recol,rep(0,length(recol)))
#dev.off()
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

#Number of objects vs. number of extinctions at the NEXT time step
R"par(mfrow=c(1,2))"
tb=1;
start=500;
R"""
plot($(obrich[start:tmax-tb]./sprich[start:tmax-tb]),$(sprich[tb+start:tmax]),pch=16,xlab='Objects per species',ylab='Number of species')
plot($(obrich[start:tmax-tb]./sprich[start:tmax-tb]),$(ext_sec[tb+start:tmax].+ext_prim[tb+start:tmax]),pch=16,xlab='Objects per species',ylab='Number of extinctions')
"""

#Ratio of objects to species vs. Richness & extincitons
R"par(mfrow=c(1,2))"
tb=1;
start=300;
rdsample = sample(collect(start:tmax-1),500,replace=false)
R"""
plot($(obrich[rdsample]./sprich[rdsample]),$(sprich[rdsample+1]),pch=16,xlab='Objects per species',ylab='Number of species')
plot($(obrich[rdsample]./sprich[rdsample]),$(ext_sec[rdsample+1].+ext_prim[rdsample+1]+1),log='y',pch=16,xlab='Objects per species',ylab='Number of extinctions')
"""

#Without Zero-inflated
R"par(mfrow=c(1,2))"
tb=1;
start=1;
n=1000;
rdsample = sample(collect(start:tmax-1),n,replace=false)
R"""
obj_ratio <- $(obrich[rdsample]./sprich[rdsample]);
sprich <- $(sprich[rdsample+1]);
ext_ratio <- $(((ext_sec[rdsample+1].+ext_prim[rdsample+1])./sprich[rdsample]));
nonzeros <- which(ext_ratio!=0);
n = length(nonzeros);
lm_rich = lm(sprich[nonzeros] ~ obj_ratio[nonzeros]);
lm_ext = lm(ext_ratio[nonzeros] ~ obj_ratio[nonzeros]);
plot(obj_ratio[nonzeros],sprich[nonzeros],pch=16,xlab='Objects per species',ylab='Number of species')
abline(lm_rich)
plot(obj_ratio[nonzeros],ext_ratio[nonzeros],pch=16,xlab='Objects per species',ylab='Extinctions per species')
abline(lm_ext)
"""



#How many cascades > a given size?
twindow = 1;
maxloss = 100;
cumext = zeros(maxloss);
for i=1:maxloss
  steps = Int(floor(tmax/twindow));
  tloss = zeros(steps);
  for t=1:steps-1
    tvec = sprich[Int(t*twindow):Int(t*twindow+twindow)]
    loss = first(tvec) - last(tvec);
    if loss < 0
      tloss[t] = 0;
    else
      tloss[t] = loss;
    end
  end
  cumext[i] = length(find(x->x>i-1,tloss));
end

R"""
plot($cumext,xlab='Extinction size',ylab='Number of extinctions',log='xy',pch=16)
"""

#
# c=sort(CID[tmax]);
# dups = Array{Int64}(0);
# for i=2:length(c)
#   if c[i]-c[i-1]==0
#     push!(dups,c[i])
#   end
# end

#Visualize int_m
int_v = Array{Int64}(length(int_m[1,:]),length(int_m[1,:]));
int_v[find(x->x=='a',int_m)]=1;
int_v[find(x->x=='n',int_m)]=2;
int_v[find(x->x=='i',int_m)]=3;
int_v[find(x->x=='m',int_m)]=4;
R"""
library(igraph)
library(plotrix)
library(RColorBrewer)
#Visualize matrix:
num_play = length($(int_v[1,:]));
plot_matrix<-function($int_v, num_play){
xx=matrix(as.numeric(as.factor($int_v)),c(num_play,num_play))
par(mar=c(1,1,1,4))
int_types=levels(as.factor($int_v))
color2D.matplot(xx,extremes=c(1:length(int_types)), border="white", axes=FALSE, xlab="", ylab="",main="")
legend(x=num_play+1,y=num_play,legend=int_types,pch=22,pt.bg=c(1:4),xpd=TRUE, bty="n")
}
plot_matrix($int_v, num_play)
"""

R"""
plot($conn,type='l',log='y',ylab='Trophic connectance')
lines($connind,lty=2)
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