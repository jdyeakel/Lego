loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

S = 400;

tmax = 1000;
tswitch = 50;

# S = 400;
probs = [
p_n=0.003,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
lambda = 0.5;

a_thresh = 0;
n_thresh = 0.2;
extmid = 0.5; #Similarity at which pr(ext) = 0.5
steep = 2; #higher is steeper

extinctions = [ones(Bool,tswitch);ones(Bool,tmax-tswitch)];
colonizations = [ones(Bool,tswitch);ones(Bool,tmax-tswitch)];

MaxN = convert(Int64,floor(S + S*lambda));

sprich = Array{Int64}(tmax);
int_m = Array{Char}();
tind_m = Array{Int64}();
prim_ext = Array{Int64}(tmax);
sec_ext = Array{Int64}(tmax);
status = Array{Int64}(tmax);
lpot_col = Array{Int64}(tmax);
CID = Array{Bool}(tmax,MaxN)*false;

maxsize = 0; tictoc=0;
#This will rerun the assembly process if the community does not assemble >10 species
while maxsize < 10
    tictoc=tictoc+1;
    println("Try ",tictoc)    
    # int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    sprich,
    cid,
    lpot_col,
    status,
    prim_ext,
    sec_ext,
    CID = assembly_trim(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
        a_thresh,n_thresh,extmid,steep,colonizations,extinctions,tmax,S,MaxN);
    
    maxsize = maximum(sum(CID,2));
    
end

R"""
par(mfrow=c(1,2))
plot($(sum(CID,2))-$sprich,type='l',lty=3,ylim=c(0,max(c($(sum(CID,2))-$sprich,$sprich))))
lines($sprich)
points($tswitch,0,pch=16,col='red')
"""

tstep = findmax(sum(CID,2))[2];
a_b,n_b,i_b,m_b,n_b0,sp_v,int_id = preamble_defs(int_m);
cid = find(isodd,CID[tstep,:]);
deg,troph = structure(S,cid,sp_v,tind_m);
spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);
#Degree distribution
# degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
adjmatrix = tind_m[[1;spcid_ind],[1;spcid_ind]];
#Visualize graph
R"""
library(igraph)
library(RColorBrewer)
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(adjmatrix'));
basal_pos <- 1
trophic = as.numeric($([0;troph[1:size(adjmatrix)[1]-1]]));
coords <- cbind(runif(vcount(fw_g)),trophic);
coords[basal_pos,1] <- 0.5
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#80808050',main=ecount(fw_g)/$(size(adjmatrix)[1])^2,vertex.label=NA,vertex.frame.color=NA, vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1)))
"""


R"plot($(sort(deg[spcid_ind],rev=true)))"


#Image the interaction matrix

#Establish community template

S = 50;
ppweight = 1/4;
# S = 400;
probs = [
# p_n=0.04,
# p_a=0.02,
# p_m=0.01,
p_n=0.05,
p_a=0.05
# p_m=0.01,
# p_i= 1 - sum([p_n,p_a]) #Ignore with 1 - pr(sum(other))
]
#expected objects per species
lambda = 0.5;

# @time int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = build_template_species(S,probs,ppweight);
# @time int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv2(S,lambda,probs,ppweight);
@time int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);


a_b,
n_b,
i_b,
m_b,
n_b0,
sp_v,
int_id = preamble_defs(int_m);

#Reorganize to clump objects
objects = deleteat!(find(x->x=='i',diag(int_m)),1);
objectssort = objects[sortperm(vec(sum(m_b[objects,objects],1)),rev=true)];
# objectssort2 = objectssort[sortperm(vec(sum(m_b[objectssort,objectssort].+a_b[objectssort,objectssort].+n_b[objectssort,objectssort],2)),rev=false)];
species = setdiff(collect(1:length(diag(int_m))),objects);
speciessort = species[sortperm(vec(sum(a_b[species,species].+n_b[species,species],1)),rev=true)];
speciessort2 = speciessort[sortperm(vec(sum(a_b[speciessort,speciessort].+n_b[speciessort,speciessort],2)),rev=false)];
int_msort = int_m[[speciessort2;objectssort],[speciessort2;objectssort]];


int_v = Array{Int64}(length(int_msort[1,:]),length(int_msort[1,:]));
int_v[find(x->x=='a',int_msort)]=1;
int_v[find(x->x=='n',int_msort)]=2;
int_v[find(x->x=='i',int_msort)]=3;
int_v[find(x->x=='m',int_msort)]=4;

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/matrixtest.pdf");
R"""
library(igraph)
library(plotrix)
library(RColorBrewer)
pal=brewer.pal(5,'Set1')
pal[4] = pal[3]
pal[3] = '#ffffff'
num_play = length($(int_v[1,:]))
xx=matrix(as.numeric(as.factor($int_v)),c(num_play,num_play))
xx2=xx;
xx2[which(xx==1)] = pal[1];
xx2[which(xx==2)] = pal[2];
xx2[which(xx==3)] = pal[3];
xx2[which(xx==4)] = pal[4];
#shade made objects
darken <- function(color, factor=1.2){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}
objects = which(diag(xx)==3)[-1];
for (i in 1:length(objects)) {
    for (j in 1:num_play) {
        col = xx2[objects[i],j];
        xx2[objects[i],j] = darken(col);
        if (length(intersect(objects,j)) == 0) {
            col = xx2[j,objects[i]];
            xx2[j,objects[i]] = darken(col);
            }
        }
    }
# pdf($namespace,height=5,width=6)
par(mar=c(1,1,1,4))
int_types=c('a','n','i','m')
color2D.matplot(xx,extremes=c(1:length(int_types)), border='white', axes=FALSE, xlab='', ylab='',main='',cellcolors=xx2)
legend(x=num_play+1,y=num_play,legend=int_types,pch=22,pt.bg=pal,xpd=TRUE, bty='n')
text(x=rep(-0.8,length(objects)),y=num_play-objects+0.5,labels='o', xpd=TRUE,cex=0.6)
text(x=objects-0.5,y=rep(num_play+0.8,length(objects)),labels='o', xpd=TRUE,cex=0.6)
text(x=-0.8,y=num_play-0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
text(x=0.8,y=num_play+0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
# dev.off()
"""
