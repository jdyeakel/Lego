using Distributions
using Gadfly
using PyCall
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/anime_sims_func.jl")





#How long does a community take to fill
rep=100;
tend = Array{Int64}(rep);

num_play = 25;
init_probs = [
p_n=0.5/num_play,
p_a=0.01,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);

tmax = num_play;
a_thresh=0.0;
n_thresh=0.2;
CA = (Array{Float64,1})[];
CAI = (Array{Float64,1})[];
CID = (Array{Int64,1})[];
for i=1:rep
  tout = 0;
  c_t = (Array{Char,1})[];
  c_tp = (Array{Int64,1})[];
  c_tind = (Array{Int64,1})[];
  cid = (Array{Int64,1})[];
  while tout==0
    tout, c_t, c_tp, c_tind, cid = anime_sims_func(int_m,tp_m,tind_m,a_thresh,n_thresh,tmax);
  end
  Cassem = Array{Float64}(tout);
  Cassemind = Array{Float64}(tout);
  for t=1:tout
    S = length(find(x->x=='n',diag(c_t[:,:,t])));
    Cassem[t] = (sum(c_tp[:,:,t])/2)/(S^2);
    Cassemind[t] = (sum(c_tind[:,:,t])/2)/(S^2);
  end
  # plot(x=collect(1:tout),y=Cassem,Scale.y_log10,Geom.point,Theme(default_color=colorant"black",default_point_size=1pt, highlight_width = 0pt))
  # plot(x=collect(1:tout),y=Cassemind,Scale.y_log10,Geom.point,Theme(default_color=colorant"black",default_point_size=1pt, highlight_width = 0pt))
  push!(CA,copy(Cassem));
  push!(CAI,copy(Cassemind));
  push!(CID,copy(cid));

  tend[i] = copy(tout);
  println(" T= ",tout)
end
vcat(CID[1:rep])


csum = (Array{Int64,1})[];
for i=1:rep
  push!(csum,cumsum(CID[i]));
end

#Visualize the assembly process
plot(
[layer(y=csum[j],x=collect(1:length(csum[j])), Geom.line, Theme(default_color=colorant"black")) for j in 1:rep]...,
Guide.xlabel("Steps"),Guide.ylabel("Species ID"))

#Connectance of each community trajectory
plot(
[layer(x=collect(1:length(CA[j])),y=CA[j], Geom.line, Theme(default_color=colorant"black")) for j in 1:rep]...,Scale.y_log10,Scale.x_log10,
Guide.xlabel("Steps"),Guide.ylabel("Species ID"))

#Difference
plot(x=collect(1:length(CID[1])),y=cumsum(CID[1])-cumsum(CID[2]),Geom.line)
