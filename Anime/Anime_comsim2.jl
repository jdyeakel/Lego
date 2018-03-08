loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

#Establish community template
S = 400;
# S = 400;
probs = [
# p_n=0.04,
# p_a=0.01,
# p_m=0.04,
p_n=0.004,
p_a=0.01,
p_m=0.002,
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

a_thresh = 0.0;
n_thresh = 0.2;
trophicload = 2;
tmax = 5000;
extinctions = true;

@time cid,
rich,
sprich,
turnover,
prim_ext,
sec_ext = assembly(a_thresh,n_thresh,trophicload,extinctions,tmax);

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
plot($degrees,$tl_ind,log='x',col=pal,pch=16,xlab="Degrees",ylab="")
"""
