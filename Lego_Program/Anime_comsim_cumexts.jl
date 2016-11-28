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

nthreshvec = [0.0,0.1,0.2];
CUMEXT = Array(Array{Int64},length(nthreshvec));
reps=10;

for n=1:length(nthreshvec)

  num_play=500;

  tmax = 1000;
  maxloss = 100;

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

  int_mv, sprich, rich, conn, comgen, ext_prim, ext_sec = repsimint(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight);


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

#
# pc = [0.8,0.6,0.4]
# R"""
# library(RColorBrewer)
# cols = brewer.pal(3,"Set1");
# plot($(CUMEXT[1][:,1]),xlab='Extinction size',ylab='Number of extinctions',log='xy',pch=16,col=cols[1],cex=$(pc[1]))
# """
# for n=1:length(nthreshvec)
#   for r=1:reps
#     R"""
#     points(jitter($(CUMEXT[n][:,r])),xlab='Extinction size',ylab='Number of extinctions',pch=16,col=cols[$n],cex=$(pc[n]))
#     """
#   end
# end
#
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/"

R"""
pdf(paste($namespace,'fig_cumulative_ext.pdf',sep=''),width=8,height=6)
library(RColorBrewer)
cols = brewer.pal(3,"Set1");
mext <- apply($(CUMEXT[1]),1,mean);
sdext <- apply($(CUMEXT[1]),1,sd);
xnums <- seq(1,100);
plot(mext,xlab='Extinction size',ylab='Number of extinctions',log='xy',type='l',col=cols[1],ylim=c(0.1,1000))
lines(mext+sdext,lty=3,col=cols[1])
lines(mext-sdext,lty=3,col=cols[1])
#polygon(x=c(xnums,rev(xnums)),y=c(mext-sdext,rev(mext+sdext)),col=paste(cols[1],50,sep=""))
"""
for n=2:length(nthreshvec)
  R"""
  mext <- apply($(CUMEXT[n]),1,mean);
  lines(mext,xlab='Extinction size',ylab='Number of extinctions',col=cols[$n])
  lines(mext+sdext,lty=3,col=cols[$n])
  lines(mext-sdext,lty=3,col=cols[$n])
  """
end
R"dev.off()"
