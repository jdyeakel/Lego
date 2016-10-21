using Distributions
using Gadfly
using PyCall



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
  int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);
  cid, c_m, crev_m, com_tp, com_tind = initiate_comm_func(int_m,tp_m,tind_m);
  for t=1:tmax
    status = "open";
    #Colonize with some probability
    rcol = rand();
    if rcol < rate_col && status == "open"
      status,cid,c_m,crev_m,com_tp,com_tind = colonize_func(int_m,tp_m,tind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
    end
    #Always run extinction code because probabilities are assessed within
    status,cid,spcid,c_m,crev_m,com_tp,com_tind = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind);
    S = length(spcid);
    # length(unique(cid))-length(cid)
    conn[t,r] = (sum(com_tp)/2)/(S^2);
    rich[t,r] = length(cid);
  end
end

#Visualize the assembly process
richplot = plot(
[layer(y=rich[:,j],x=collect(1:tmax), Geom.line, Theme(default_color=colorant"gray")) for j in 1:reps]...,
Guide.xlabel("Time"),Guide.ylabel("Richness"),Scale.x_log10);

draw(PDF("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_fullassemb.pdf", 5inch, 4inch), assembplot)





#THIS WORKS! 10/21/16
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
num_play = 500;
init_probs = [
p_n=1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
a = SharedArray(Int64,10)
@sync @parallel for i=1:10
  int_m, sp_m, t_m, tp_m, tind_m = build_template_degrees(num_play,init_probs);
  a[i] = length(find(x->x=='n',int_m));
end


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
num_play = 500;
init_probs = [
p_n=1/num_play,
p_a=0.01,
p_m=0.1/num_play,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]
a = @spawn build_template_degrees(num_play,init_probs);
b = @spawn build_template_degrees(num_play,init_probs);
[fetch(a), fetch(b)]


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/count_heads.jl")

a = @spawn count_heads(100000000)
b = @spawn count_heads(100000000)
fetch(a)+fetch(b)
