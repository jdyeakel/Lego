loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

#Establish community template
S = 15;
# S = 400;
probs = [
p_n=0.04,
p_a=0.01,
p_m=0.04,
# p_n=0.004,
# p_a=0.01,
# p_m=0.002,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]

ppweight = 1/4;
sim=false;
par=false;
calcpotcol = false;
@time int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);

a_b,
n_b,
i_b,
m_b,
n_b0,
sp_v,
int_id = preamble_defs(int_m);


###################################
# COMMUNITY SIMULATION THROUGH TIME
###################################
#Establish colonization and extinction rates

cid = Array{Int64}(0);
a_thresh = 0.0;
n_thresh = 0.2;
trophicload = 2;
tmax = 2000;
extinctions = true;

rich = Array{Int64}(tmax);
sprich = Array{Int64}(tmax);
obrich = Array{Int64}(tmax);

@time for t = 1:tmax
  #Print every 1000 timesteps
  if mod(t,1000)==0
    println(string("t=",t))
  end

  cid,
  lpot_col,
  status,
  num_ext1,
  num_ext2 = colext(int_m,cid,a_thresh,n_thresh,extinctions,trophicload);
  
  spcid = intersect(sp_v,cid);
  spcid_ind = indexin(spcid,[1;sp_v]);
  
  rich[t] = length(cid);
  sprich[t] = length(spcid);
  obrich[t] = rich[t] - sprich[t];

end

spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);

degrees,tl_ind,conn = structure();

R"par(mfrow=c(1,2))"
R"""
plot($(collect(1:tmax)),$rich,type='l',lty=2,log="x",xlab="Time",ylab="Richness")
lines($(collect(1:tmax)),$sprich)
"""

R"""
library(RColorBrewer)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($degrees))
plot($degrees,$tl_ind,log='x',col=pal,pch=16)
"""
