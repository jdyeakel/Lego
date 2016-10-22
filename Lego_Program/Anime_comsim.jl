using Distributions
using Gadfly
using PyCall
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")

#Establish community template
num_play = 500;
init_probs = [
p_n=1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
sim=true;
int_m, sp_m, t_m, tp_m, tind_m, simvalue = build_template_degrees(num_play,init_probs,sim);

#Establish colonization and extinction rates
rate_col = 0.2;
rate_ext = 0.1;

#Establish thresholds
a_thresh = 0.2;
n_thresh = 0.4;

rep = 1000;
CID = (Array{Int64,1})[];
rich = Array{Int64}(rep);
for r = 1:rep
  #The add-until-full simulation
  #Creating a new int_m each time
  sim=true;
  int_m, sp_m, t_m, tp_m, tind_m, simvalue = build_template_degrees(num_play,init_probs,sim);
  cid, c_m, crev_m, com_tp, com_tind = initiate_comm_func(int_m);
  status = "open";
  while status == "open"
    status,cid,c_m,crev_m,com_tp,com_tind = colonize_func(a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
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
num_play = 500;
init_probs = [
p_n=1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
#int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);

#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0;
n_thresh = 0.2;
tmax = 1000;
CID = (Array{Int64,1})[];
rich = Array{Int64}(tmax);
conn = Array{Float64}(tmax);
sim=true;
int_m, sp_m, t_m, tp_m, tind_m, simvalue = build_template_degrees(num_play,init_probs,sim);
cid, c_m, crev_m, com_tp, com_tind = initiate_comm_func(int_m,tp_m,tind_m);
@time for t = 1:tmax
  #The add-until-full simulation
  #Creating a new int_m each time
  status = "open";
  #Colonize with some probability
  rcol = rand();
  if rcol < rate_col && status == "open"
    status,cid,c_m,crev_m,com_tp,com_tind = colonize_func(int_m,tp_m,tind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
  end
  #Always run extinction code because probabilities are assessed within
  status,cid,spcid,c_m,crev_m,com_tp,com_tind = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,simvalue);
  S = length(spcid);
  # length(unique(cid))-length(cid)
  conn[t] = (sum(com_tp)/2)/(S^2);
  rich[t] = length(cid);
  println("Richness = ",rich[t])
  push!(CID,copy(cid));
end
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
tmax = 1000;
reps=50;

#Shared variables
sprich = SharedArray{Int64}(tmax,reps);
rich = SharedArray{Int64}(tmax,reps);
conn = SharedArray{Float64}(tmax,reps);

@sync @parallel for r=1:reps
  #Establish community template
  num_play = 500;
  init_probs = [
  p_n=1/num_play,
  p_a=0.01,
  p_m=0.1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]
  sim=true;
  int_m, sp_m, t_m, tp_m, tind_m, simvalue = build_template_degrees(num_play,init_probs,sim);
  cid, c_m, crev_m, com_tp, com_tind = initiate_comm_func(int_m,tp_m,tind_m);
  for t=1:tmax
    status = "open";
    #Colonize with some probability
    rcol = rand();
    if rcol < rate_col && status == "open"
      status,cid,c_m,crev_m,com_tp,com_tind = colonize_func(int_m,tp_m,tind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
    end
    #Always run extinction code because probabilities are assessed within
    status,cid,spcid,c_m,crev_m,com_tp,com_tind = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,simvalue);
    sprich[t,r] = length(spcid);
    # length(unique(cid))-length(cid)
    conn[t,r] = (sum(com_tp)/2)/(sprich[t,r]^2);
    rich[t,r] = length(cid);
  end
end

#Visualize richness over time
richplot = plot(
[layer(y=sprich[:,j],x=collect(1:tmax), Geom.line, Theme(default_color=colorant"gray")) for j in 1:reps]...,
Guide.xlabel("Time"),Guide.ylabel("Richness"),Scale.y_log10);
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_diversity.pdf", 5inch, 4inch), richplot)

#Connectance Plot
connplot = plot(
[layer(y=conn[:,j],x=collect(1:tmax), Geom.line, Theme(default_color=colorant"gray")) for j in 1:reps]...,
Guide.xlabel("Time"),Guide.ylabel("Richness"),Scale.x_log10,Scale.y_log10);
draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_connectance.pdf", 5inch, 4inch), connplot)
