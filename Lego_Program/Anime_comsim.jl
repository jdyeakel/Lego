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
init_probs = [
p_n=0.1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
sim=true;
ppweight=1/3;
int_m, sp_m, t_m, tp_m, tind_m, mp_m, simvalue = build_template_degrees(num_play,init_probs,ppweight,sim);

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
  int_m, sp_m, t_m, tp_m, tind_m, mp_m, simvalue = build_template_degrees(num_play,init_probs,ppweight,sim);
  cid, c_m, crev_m, com_tp, com_tind = initiate_comm_func(int_m,tp_m,tind_m);
  status = "open";
  while status == "open"
    status,cid,c_m,crev_m,com_tp,com_tind = colonize_func(int_m,tp_m,tind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
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


#Establish community template
num_play = 20;
init_probs = [
p_n=0.1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
#int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);

#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0.0;
n_thresh = 0.2;
tmax = 1000;
CID = (Array{Int64,1})[];
rich = Array{Int64}(tmax);
conn = Array{Float64}(tmax);
ext_prim = Array{Int64}(tmax);
ext_sec = Array{Int64}(tmax);
comgen =zeros(Int64,tmax,num_play);
ppweight = 1/3;
sim=true;
int_m, sp_m, t_m, tp_m, tind_m, mp_m, simvalue = build_template_degrees(num_play,init_probs, ppweight, sim);
cid, c_m, crev_m, com_tp, com_tind = initiate_comm_func(int_m,tp_m,tind_m);
@time for t = 1:tmax
  #The add-until-full simulation
  #Creating a new int_m each time
  #status = "open"; #I don't think we need this now
  #Colonize with some probability
  rcol = rand();
  if rcol < rate_col #&& status == "open"
    status,cid,c_m,crev_m,com_tp,com_tind = colonize_func(int_m,tp_m,tind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
  end
  #Always run extinction code because probabilities are assessed within
  status,cid,spcid,c_m,crev_m,com_tp,com_tind,extinctions = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,simvalue);
	#Save primary and secondary extinction information
	ext_prim[t] = extinctions[1];
	ext_sec[t] = extinctions[2];
  S = length(spcid);
  # length(unique(cid))-length(cid)
  conn[t] = (sum(com_tp))/(S^2);
  rich[t] = length(cid);
  println("Richness = ",rich[t])
  push!(CID,copy(cid));
  comgen[t,cid] = 1;
end
writedlm("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen.csv",comgen);


csum = (Array{Int64,1})[];
for i=1:rep
  push!(csum,cumsum(CID[i]));
end
#Diversity through time
divplot = plot(y=rich,x=collect(1:tmax), Geom.line, Scale.x_log10, Theme(default_color=colorant"black"),
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
reps=200;

num_play = 500;
#Shared variables
sprich = SharedArray{Int64}(tmax,reps);
rich = SharedArray{Int64}(tmax,reps);
conn = SharedArray{Float64}(tmax,reps);
comgen = SharedArray{Int64}(reps,num_play,tmax);
ext_prim = SharedArray{Int64}(tmax,reps);
ext_sec = SharedArray{Int64}(tmax,reps);

init_probs = [
p_n=1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
ppweight = 1/3;
sim=true;
int_m, sp_m, t_m, tp_m, tind_m, mp_m, simvalue = build_template_degrees(num_play,init_probs,ppweight,sim);

@sync @parallel for r=1:reps
  #Establish community template
  cid, c_m, crev_m, com_tp, com_tind = initiate_comm_func(int_m,tp_m,tind_m);
  for t=1:tmax
    status = "open";
    #Colonize with some probability
    rcol = rand();
    if rcol < rate_col && status == "open"
      status,cid,c_m,crev_m,com_tp,com_tind = colonize_func(int_m,tp_m,tind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
    end
    #Always run extinction code because probabilities are assessed within
    status,cid,spcid,c_m,crev_m,com_tp,com_tind,extinctions = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,simvalue);
    #Save primary and secondary extinction information
  	ext_prim[t,r] = extinctions[1];
  	ext_sec[t,r] = extinctions[2];
    sprich[t,r] = length(spcid);
    # length(unique(cid))-length(cid)
    conn[t,r] = (sum(com_tp)/2)/(sprich[t,r]^2);
    rich[t,r] = length(cid);
    comgen[r,cid,t] = 1;
  end
end

#Similarity posthoc analysis
jacc = Array{Float64}(reps,reps,tmax);
meansim = Array{Float64}(reps,tmax);
maxsim = Array{Float64}(reps,tmax);
for t=1:tmax
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

#Visualize richness over time
tmax = 500
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
