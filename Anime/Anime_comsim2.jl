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
tmax = 100;
extinctions = true;


@time cid,
rich,
sprich,
turnover,
mres_overlap,
res_overlap_dist,
conn,
conn_ind,
prim_ext,
sec_ext,
status,
lpot_col,
avgdegree,
degrees,
trophic = assembly(
    int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
    a_thresh,n_thresh,extinctions,tmax,S);

spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);
degrees,tl_ind = structure(S,cid,sp_v,tind_m);
#Null degrees,tl_ind
cid_null = collect(1:size(int_m)[1]);
degrees_null,tl_ind_null = structure(S,cid_null,sp_v,tind_m);

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/res_overlap.pdf");
R"""
pdf($namespace,height=5,width=6)
plot($res_overlap,log='x',xlab="Time",ylab="Resource overlap",type='l')
dev.off()
"""
namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/richness_trophic.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=10)
par(mfrow=c(1,2))
plot($(collect(1:tmax)),$rich,type='l',lty=2,log="x",xlab="Time",ylab="Richness")
lines($(collect(1:tmax)),$sprich)
pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($degrees))
plot($degrees_null,$tl_ind_null,log='xy',col="808080",pch=16,xlab="Degrees",ylab="Trophic level",ylim=c(1,max(cbind($tl_ind,$tl_ind_null))))
points($degrees,$tl_ind,col=pal,pch=16)
dev.off()
"""
