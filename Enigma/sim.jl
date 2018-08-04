loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly2.jl")


S = 400;

maxits = 1000;

# S = 400;
probs = [
p_n=0.001,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
lambda = 0.5;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));

int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

a_b,
n_b,
i_b,
m_b,
n_b0,
sp_v,
int_id = preamble_defs(int_m);

@time sprich,rich,clock,CID = assembly(
    int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
    athresh,nthresh,maxits);

# @time sprich,rich,clock = assembly2(
#     int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
#     athresh,nthresh,tmax);

namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/sprich_web.pdf"
R"""
pdf($namespace,width=10,height=5)
par(mfrow=c(1,2))
plot($clock,$sprich,type='l',xlab='Time',ylab='Sp/Ob richness',ylim=c(0,max($([sprich;rich.-sprich]))))
lines($clock,$(rich .- sprich),col='gray')
#dev.off()
"""

tstep = maxits;
cid = find(isodd,CID[:,tstep]);
deg,troph = structure(S,cid,sp_v,tind_m);
spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);
#Degree distribution
# degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
adjmatrix = tind_m[[1;spcid_ind],[1;spcid_ind]];
indmatrix = adjmatrix .- tp_m[[1;spcid_ind],[1;spcid_ind]];
dirmatrix = tp_m[[1;spcid_ind],[1;spcid_ind]];

#Visualize graph
namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/webnet.pdf"
R"""
library(igraph)
library(RColorBrewer)
#pdf($namespace,width=6,height=5)
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(adjmatrix'));
basal_pos <- 1
trophic = as.numeric($([0;troph[1:size(adjmatrix)[1]-1]]));
keepnodes = c(1,which(trophic>0.1))"""; @rget keepnodes; keepnodes = Int64.(keepnodes);
R"""
trophic2 = trophic[keepnodes]
coords <- cbind(runif(length(keepnodes)),trophic2);
coords[basal_pos,1] <- 0.5
fw_g = graph.adjacency($(adjmatrix[keepnodes,keepnodes]'))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1))) 
#main=ecount(fw_g)/$(size(adjmatrix)[1])^2,
fw_ind <- graph.adjacency($(indmatrix[keepnodes,keepnodes]'));
#plot(fw_ind,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='red',vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)),add=TRUE)
dev.off()
"""
