using Distributions
using Gadfly
using PyCall
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")

#Establish community template
num_play = 100;
init_probs = [
p_n=0.5/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);

#Establish colonization and extinction rates
rate_col = 0.2;
rate_ext = 0.1;

#Establish thresholds
a_thresh = 0.2;
n_thresh = 0.4;

rep = 100;
CID = (Array{Int64,1})[];
for r = 1:rep
  #The add-until-full simulation
  cid, c_m, crev_m, com_sparse, com_tp, com_tind = initiate_comm_func(int_m);
  status = "open";
  while status == "open"
    status, cid, c_m, crev_m, com_sparse, com_tp, com_tind = colonize_func(a_thresh,n_thresh,cid,c_m,crev_m,com_sparse,com_tp,com_tind);
  end
  rich = length(cid);
  println("Richness = ",rich)
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
