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
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsimint.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonize_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")

nthreshvec = [0.05,0.1,0.2];


for n=1:length(nthreshvec)
  
  num_play=500;
  reps=20;
  tmax = 5000;
  maxloss = 100;
  CUMEXT = Array(Array{Int64},length(nthreshvec));
  ppweight = 1/4;

  a_thresh=0.0;
  n_thresh = nthreshvec[n];
  trophicload = 2;
  rate_col=1;
  
  #Establish community template
  probs = [
  p_n=0.01,
  p_a=0.01,
  p_m= 0.002, #needvec[i]/num_play, #1/num_play,
  p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
  ]
  
  int_m, sprich, rich, conn, comgen, ext_prim, ext_sec = repsimint(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);

  
  #How many cascades > a given size?
  cumext = SharedArray{Int64}(maxloss,reps);
  @sync @parallel for r=1:reps
    twindow = 1;
    sprichr = sprich[:,r];
    for i=1:maxloss
      steps = Int(floor(tmax/twindow));
      tloss = zeros(steps);
      for t=1:steps-1
        tvec = sprichr[Int(t*twindow):Int(t*twindow+twindow)]
        loss = first(tvec) - last(tvec);
        if loss < 0
          tloss[t] = 0;
        else
          tloss[t] = loss;
        end
      end
      cumext[i,r] = length(find(x->x>i-1,tloss));
    end
  end
  CUMEXT[n] = cumext;
end



R"""
library(RColorBrewer)
cols = brewer.pal(3,"Set1");
plot($(CUMEXT[1][:,1]),xlab='Extinction size',ylab='Number of extinctions',log='xy',pch=16,col=cols[1],cex=0.5)
"""
for n=1:length(nthreshvec)
  for r=1:reps
    R"""
    points($(CUMEXT[n][:,r]),xlab='Extinction size',ylab='Number of extinctions',log='xy',pch=16,col=cols[$n],cex=0.5)
    points($(CUMEXT[n][:,r]),xlab='Extinction size',ylab='Number of extinctions',log='xy',pch=16,col=cols[$n],cex=0.5)
    """
  end
end
