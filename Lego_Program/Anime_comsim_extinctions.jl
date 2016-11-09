
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

needvec = collect(0.5:0.5:5);
ln = length(needvec);

RvsE = Array(Array{Float64},ln);
rich_ss = zeros(Float64,ln);
#Loop over need value
for i=1:ln
  
  num_play = 500;
  #Shared variables
  sprich = SharedArray{Int64}(tmax,reps);
  rich = SharedArray{Int64}(tmax,reps);
  conn = SharedArray{Float64}(tmax,reps);
  comgen = SharedArray{Int64}(reps,num_play,tmax);
  ext_prim = SharedArray{Int64}(tmax,reps);
  ext_sec = SharedArray{Int64}(tmax,reps);
  
  #Establish community template
  probs = [
  p_n=needvec[i]/num_play,
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
  
  #Number of extinctions for a community of a given species richness
  mext = zeros(maximum(sprich),2);
  for j=1:maximum(sprich)
    numsp=find(x->x==j,sprich);
    lnumsp = length(numsp);
    randsamp = sample(collect(1:lnumsp),100,replace=true);
    mext[j,1] = j;
    mext[j,2] = mean(ext_sec[numsp][randsamp]+ext_prim[numsp][randsamp]);
  end
  RvsE[i] = mext;
  rich_ss[i] = mean(sprich[(tmax-100):(tmax),:]);
end


#Plot Extinctions vs. Richness
R"""
library(RColorBrewer);
cols <- brewer.pal(10,'Spectral');
x=$(RvsE[1][:,1]);
y=$(RvsE[1][:,2]);
plot(x,y,type='l',ylim=c(0,15),col=cols[1],
  xlab="Species richness",ylab="Number of extinctions");
points($(rich_ss[1]),10,pch=16,col=cols[1]);
segments($(rich_ss[1]),10,$(rich_ss[1]),0,cols[1]);
"""
for i=2:10
  x=(RvsE[i][:,1]);
  y=(RvsE[i][:,2]);
  richss=rich_ss[i];
  R"lines($x,$y,col=cols[$i]);
  points($richss,10,pch=16,col=cols[$i])
  "
end






richextplot = plot(
[layer(x=RvsE[j][:,1],y=RvsE[j][:,2], Geom.line, Theme(default_color=colorant"black",default_point_size=2pt, highlight_width = 0pt)) for j in 1:ln]...,
Guide.xlabel("Species richness"),Guide.ylabel("Mean total extinctions"));
