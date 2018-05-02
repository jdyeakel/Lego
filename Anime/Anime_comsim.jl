loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

# #Establish community template
# S = 20;
# lambda = 0.25;
# ppweight = 1/4;
# # S = 400;
# probs = [
# p_n=0.04,
# p_a=0.08,
# p_m=0.01,
# # p_n=0.004,
# # p_a=0.01,
# # p_m=0.002,
# p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
# ]
# 

S = 400;
ppweight = 1/4;
# S = 400;
probs = [
p_n=0.01,
p_a=0.01
]
#expected objects per species
lambda = 1;

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


###################################
# COMMUNITY SIMULATION THROUGH TIME
###################################
#Establish colonization and extinction rates
a_thresh = 0.0;
n_thresh = 0.2;
tmax = 1000;
tswitch = 500;
extinctions = [ones(Bool,tswitch);ones(Bool,tswitch)];
colonizations = [ones(Bool,tswitch);ones(Bool,tswitch)];

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
    a_thresh,n_thresh,colonizations,extinctions,tmax,S);

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures2/rich_time.pdf");
R"""
pdf($namespace,height=5,width=6)
plot($sprich,type='l',ylim=c(0,max(c($(sprich),$rich-$sprich))),xlab='Time',ylab='Richness')
lines($rich-$sprich,lty=2)
points($tswitch,0,pch=16,col='red')
dev.off()
"""


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






loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

S = 400;
ppweight = 1/4;
# S = 400;
probs = [
p_n=0.0025,
p_a=0.0025
]
#expected objects per species
lambda = 0.5;

a_thresh = 0;
n_thresh = 0.2;
tmax = 2000;
tswitch = 1000;
extinctions = [ones(Bool,tswitch);ones(Bool,tmax-tswitch)];
colonizations = [ones(Bool,tswitch);zeros(Bool,tmax-tswitch)];

MaxN = convert(Int64,floor(S + S*lambda));

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

    cid,
    lpot_col,
    status,
    prim_ext,
    sec_ext,
    CID = assembly_trim(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
        a_thresh,n_thresh,colonizations,extinctions,tmax,S,MaxN);
    
    maxsize = maximum(sum(CID,2));
    
end
R"""plot($(sum(CID,2)),type='l')
points($tswitch,0,pch=16,col='red')"""












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
