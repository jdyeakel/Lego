using Distributions
using Gadfly
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/anime_sims_func.jl")





#How long does a community take to fill
rep=50;
tend = Array{Int64}(rep);

num_play = 100;
init_probs = [
p_n=0.01,
p_a=0.01,
p_m=0.001,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
a_thresh=0.0;
n_thresh=0.01;
CA = (Array{Float64,1})[];
CAI = (Array{Float64,1})[];
int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);
for i=1:rep
  tout = 0;
  Cassem = Array{Float64}(tout);
  Cassemind = Array{Float64}(tout);
  while tout==0
    tout, c_t, c_tp, c_tind = anime_sims_func(int_m,tp_m,tind_m,a_thresh,n_thresh);


    for t=1:tout
      S = length(find(x->x=='n',diag(c_t[:,:,t])));
      Cassem[t] = (sum(c_tp[:,:,t])/2)/(S^2);
      Cassemind[t] = (sum(c_tind[:,:,t])/2)/(S^2);
    end

  end
  # plot(x=collect(1:tout),y=Cassem,Scale.y_log10,Geom.point,Theme(default_color=colorant"black",default_point_size=1pt, highlight_width = 0pt))
  # plot(x=collect(1:tout),y=Cassemind,Scale.y_log10,Geom.point,Theme(default_color=colorant"black",default_point_size=1pt, highlight_width = 0pt))
  CA[i] = copy(Cassem);
  CAI[i] = copy(Cassemind);
  
  tend[i] = copy(tout);
  print("Rep= ",i)
end
plot(x=tend,Geom.histogram)
