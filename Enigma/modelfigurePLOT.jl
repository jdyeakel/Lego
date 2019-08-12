if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end


S = 200;

maxits =100;
cn = pi;
ce = sqrt(2);
cp = 1;
athresh = 0.0;
nthresh = 1.0;

SOprobs = (
p_n=0.002,
p_a=0.01
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


#expected objects per species
lambda = 0.2;
MaxN = convert(Int64,floor(S + S*lambda));

int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

a_b,n_b,i_b,m_b,n_b0,sp_v,int_id = preamble_defs(int_m);

#Build SS trophic net
@time sprich,rich,clock,CID = assembly(
    int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
    athresh,nthresh,maxits,cn,ce,cp);

#Reorganize to clump objects
# objects = deleteat!(findall(x->x=='i',diag(int_m)),1);
# objectssort = objects[sortperm(vec(sum(m_b[objects,objects],dims=1)),rev=true)];




# 
# tstep = maxits;
# cid = findall(isodd,CID[:,tstep]);
# deg,troph = structure(S,cid,sp_v,tind_m);
# spcid = intersect(sp_v,cid);
# spcid_ind = indexin(spcid,[1;sp_v]);
# #Degree distribution
# # degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
# adjmatrix = tind_m[[1;spcid_ind],[1;spcid_ind]];
# indmatrix = adjmatrix .- tp_m[[1;spcid_ind],[1;spcid_ind]];
# dirmatrix = tp_m[[1;spcid_ind],[1;spcid_ind]];
# 
# species = [1;spcid];
# objects = setdiff(cid,spcid);
# num_play = length(species) + length(objects);
# 
# # objects = deleteat!(findall(x->x=='i',diag(int_m)),1);
# objectssort = objects[sortperm(vec(sum(m_b[objects,objects],dims=1)),rev=true)];
# # species = setdiff(collect(1:length(diag(int_m))),objects);
# speciessort = species[sortperm(vec(sum(a_b[species,species].+n_b[species,species],dims=1)),rev=true)];
# speciessort2 = speciessort[sortperm(vec(sum(a_b[speciessort,speciessort].+n_b[speciessort,speciessort],dims=2)),rev=false)];
# int_msort = int_m[[speciessort2;objectssort],[speciessort2;objectssort]];
# 
# # int_msort = int_m[[species;objects],[species;objects]];
# int_v = Array{Int64}(undef,length(int_msort[1,:]),length(int_msort[1,:]));
# int_v[(LinearIndices(int_msort))[findall(x->x=='a',int_msort)]].=1;
# int_v[(LinearIndices(int_msort))[findall(x->x=='n',int_msort)]].=2;
# int_v[(LinearIndices(int_msort))[findall(x->x=='i',int_msort)]].=3;
# int_v[(LinearIndices(int_msort))[findall(x->x=='m',int_msort)]].=4;
# 
# namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/matrixtest.pdf");
# R"""
# library(igraph)
# library(plotrix)
# library(RColorBrewer)
# library(SDMTools)
# library(ellipse)
# pal=brewer.pal(5,'Set1')
# pal[4] = pal[3]
# pal[3] = '#ffffff'
# num_play = length($(int_v[1,:]))
# xx=matrix(as.numeric(as.factor($int_v)),c(num_play,num_play))
# xx2=xx;
# xx2[which(xx==1)] = pal[1];
# xx2[which(xx==2)] = pal[2];
# xx2[which(xx==3)] = pal[3];
# xx2[which(xx==4)] = pal[4];
# #shade made objects
# darken <- function(color, factor=1.2){
#     col <- col2rgb(color)
#     col <- col/factor
#     col <- rgb(t(col), maxColorValue=255)
#     col
# }
# objects = which(diag(xx)==3); #[-1]
# for (i in 1:length(objects)) {
#     for (j in 1:num_play) {
#         col = xx2[objects[i],j];
#         xx2[objects[i],j] = darken(col);
#         if (length(intersect(objects,j)) == 0) {
#             col = xx2[j,objects[i]];
#             xx2[j,objects[i]] = darken(col);
#             }
#         }
#     }
# pdf($namespace,height=5,width=10)
# par(mfrow=c(1,2))
# #pdf($namespace,width=6,height=5)
# pal <- brewer.pal(3,"Set1")
# fw_g <- graph.adjacency($(transpose(adjmatrix)));
# basal_pos <- 1
# trophic = as.numeric($([0;troph[1:size(adjmatrix)[1]-1]]));
# # trophic = as.numeric($([0;paths[keepnodes[2:length(keepnodes)]]]));
# keepnodes = c(1,which(trophic>0.9))
# """
# @rget keepnodes; keepnodes = Int64.(keepnodes);
# R"""
# #keepnodes = $keepnodes;
# trophic2 = trophic[keepnodes];
# #BUILD TROPHIC COORDS TO GRAB FROM
# maxTL = floor(max(trophic2)) + 1
# TLcoords = list()
# for (i in 1:maxTL) {
#     whole_ell = ellipse(matrix(c(0,0.0,0.0,0),2,2),c(2/(i+5),1))
#     low = which(whole_ell[,2] < 1 - i*0.2)
#     TLcoords[[i]] = whole_ell[low,]
#     TLcoords[[i]][,2] = TLcoords[[i]][,2] + i + 2 + rnorm(length(TLcoords[[i]][,2]),0,0.1)
# 
# }
# coords = matrix(0,nrow=length(trophic2),ncol=2)
# for (i in 2:length(trophic2)) {
#     trfloor = floor(trophic2[i])
#     if (trfloor < 1) {
#         trfloor = 1
#     }
#     #grab random position
#     s = sample(seq(1:length(TLcoords[[trfloor]][,2])),1)
#     coords[i,] = TLcoords[[trfloor]][s,]
# }
# # coords <- cbind(runif(length(keepnodes)),trophic2);
# coords[basal_pos,] <- c(0,-1)
# max_pos = which(trophic2 == max(trophic2))[1]
# coords[max_pos,1] <- 0
# objectcoords <- cbind(seq(-1,1,length.out=length(objects)),rep(maxTL+2,length(objects)));
# colpos = floor(trophic2)+1
# lowtrophic = which(colpos == 1)
# colpos[lowtrophic[2:length(lowtrophic)]] = 2
# nodecols = c("white",rev(colorRampPalette(brewer.pal(9,"Spectral"))(max(colpos))))
# fw_g = graph.adjacency($(transpose(adjmatrix[keepnodes,keepnodes])))
# # plot(fw_g,layout=layout_(fw_g,nicely()),vertex.size=4,arrow.size=0.25)
# plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color='#6495ED',vertex.label=NA, vertex.color=nodecols[colpos],vertex.frame.color="black",cex=2)
# legend.gradient(cbind(x =c(1.2,1.3,1.3,1.2), y =c(-.3,-.3,0.8,0.8)),cols=nodecols,limits = c(0,floor(max(trophic2)))+1,title='TL')
# #vertex.color=c(pal[1],rep(pal[2],vcount(fw_g)-1))
# #main=ecount(fw_g)/$(size(adjmatrix)[1])^2,
# fw_ind <- graph.adjacency($(transpose(indmatrix[keepnodes,keepnodes])));
# # plot(fw_ind,layout=coords,vertex.size=5,edge.arrow.size=0.25,edge.color=pal[1],vertex.label=NA, vertex.color=nodecols[colpos],vertex.frame.color="black",cex=2,add=TRUE)
# # dev.off()
# fw_obe <- graph.adjacency($(transpose(a_b[[keepnodes;objects],[keepnodes;objects]])));
# plot(fw_obe,layout=rbind(coords,objectcoords),vertex.size=5,edge.arrow.size=0.25,edge.color=pal[1],vertex.label=NA, vertex.color=c(nodecols[colpos],rep('black',length(objects))),vertex.frame.color="black",cex=2,add=TRUE)
# int_types=c('e','n','i','m')
# color2D.matplot(xx,extremes=c(1:length(int_types)), border='white', axes=FALSE, xlab='', ylab='',main='',cellcolors=xx2)
# legend(x=-15,y=num_play,legend=int_types,pch=22,pt.bg=pal,xpd=TRUE, bty='n')
# dev.off()
# """
# 











#PLOT SPECIES AND OBJECTS


tstep = maxits;
cid = findall(isodd,CID[:,tstep]);
deg,troph = structure(S,cid,sp_v,tind_m);
spcid = intersect(sp_v,cid);
spcid_ind = indexin(spcid,[1;sp_v]);
#Degree distribution
# degrees = vec(sum(tind_m[spcid_ind,spcid_ind],2));
adjmatrix = a_b[[1;cid],[1;cid]];
adjmatrix_n = n_b0[[1;cid],[1;cid]];
adjmatrix_m = m_b[[1;cid],[1;cid]];
# indmatrix = adjmatrix .- tp_m[[1;spcid_ind],[1;spcid_ind]];
# dirmatrix = tp_m[[1;spcid_ind],[1;spcid_ind]];

species = [1;spcid];
objects = setdiff(cid,spcid);
num_play = length(species) + length(objects);

# objects = deleteat!(findall(x->x=='i',diag(int_m)),1);
objectssort = objects[sortperm(vec(sum(m_b[objects,objects],dims=1)),rev=true)];
# species = setdiff(collect(1:length(diag(int_m))),objects);
speciessort = species[sortperm(vec(sum(a_b[species,species].+n_b[species,species],dims=1)),rev=true)];
speciessort2 = speciessort[sortperm(vec(sum(a_b[speciessort,speciessort].+n_b[speciessort,speciessort],dims=2)),rev=false)];
int_msort = int_m[[speciessort2;objectssort],[speciessort2;objectssort]];

# int_msort = int_m[[species;objects],[species;objects]];
int_v = Array{Int64}(undef,length(int_msort[1,:]),length(int_msort[1,:]));
int_v[(LinearIndices(int_msort))[findall(x->x=='a',int_msort)]].=1;
int_v[(LinearIndices(int_msort))[findall(x->x=='n',int_msort)]].=2;
int_v[(LinearIndices(int_msort))[findall(x->x=='i',int_msort)]].=3;
int_v[(LinearIndices(int_msort))[findall(x->x=='m',int_msort)]].=4;

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/figures/matrixtest.pdf");
R"""
library(igraph)
library(plotrix)
library(RColorBrewer)
library(SDMTools)
library(ellipse)
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
objects = which(diag(xx)==3);
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
pdf($namespace,height=5,width=10)
par(mfrow=c(1,2))
pal <- brewer.pal(3,"Set1")
fw_g <- graph.adjacency($(transpose(adjmatrix)));
basal_pos <- 1
trophic = as.numeric(c($([0;troph[1:length(species)-1]]),seq(1,8,length.out=length($objects))));
keepnodes = c(1,which(trophic>0.0))
"""
@rget keepnodes; keepnodes = Int64.(keepnodes);
R"""
trophic2 = trophic[keepnodes];
maxTL = floor(max(trophic2)) + 1
TLcoords = list()
for (i in 1:maxTL) {
    whole_ell = ellipse(matrix(c(0,0.0,0.0,0),2,2),c(2/(2*sqrt(i)+4),1))
    low = which(whole_ell[,2] < 1 - i*0.2)
    TLcoords[[i]] = whole_ell[low,]
    TLcoords[[i]][,2] = TLcoords[[i]][,2] + i + 2 + rnorm(length(TLcoords[[i]][,2]),0,0.1)
    
}
coords = matrix(0,nrow=length(trophic2),ncol=2)
for (i in 2:(length(trophic2)-length(objects)+1)) {
    trfloor = floor(trophic2[i])
    if (trfloor < 1) {
        trfloor = 1
    }
    s = sample(seq(1:length(TLcoords[[trfloor]][,2])),1)
    coords[i,] = TLcoords[[trfloor]][s,]
}
for (i in ((length(trophic2)-length(objects)+2):length(trophic2))) {
    coords[i,] = c(-1.2,trophic2[i])
}
coords[basal_pos,] <- c(0,-1)
max_pos = which(trophic2[1:length($species)] == max(trophic2[1:length($species)]))[1]
coords[max_pos,1] <- 0
coords[(length(trophic2)-length($objects)+1):length(trophic2),2] = seq(min(coords[,2]),max(coords[1:(length(trophic2)-length($objects)),2]),length.out=length($objects))
colpos = floor(trophic2)+1
lowtrophic = which(colpos == 1)
colpos[lowtrophic[2:length(lowtrophic)]] = 2
nodecols = c("white",rev(colorRampPalette(brewer.pal(9,"Spectral"))(max(colpos))))
vertexcols = nodecols[colpos];
vertexcols[((length(trophic2)-length(objects)+2):length(trophic2))] ="black";
fw_g = graph.adjacency($(transpose(adjmatrix[keepnodes,keepnodes])))
plot(fw_g,layout=coords,vertex.size=8,edge.width=2,edge.arrow.size=0.5,edge.color=paste(pal[1],'75',sep=''),vertex.label=NA, vertex.color=vertexcols,vertex.frame.color="black",cex=2)
legend.gradient(cbind(x =c(1.2,1.3,1.3,1.2), y =c(-.3,-.3,0.8,0.8)),cols=nodecols,limits = c(0,floor(max(trophic2)))+1,title='TL')
fw_needs <- graph.adjacency($((adjmatrix_n[keepnodes,keepnodes])));
plot(fw_needs,layout=coords,vertex.size=8,edge.width=2,edge.arrow.size=0.5,edge.color=paste(pal[2],'75',sep=''),vertex.label=NA, vertex.color=vertexcols,vertex.frame.color="black",cex=2,add=TRUE)
fw_makes <- graph.adjacency($((adjmatrix_m[keepnodes,keepnodes])));
plot(fw_makes,layout=coords,vertex.size=8,edge.width=2,edge.arrow.size=0.5,edge.color=paste(pal[3],'75',sep=''),vertex.label=NA, vertex.color=vertexcols,vertex.frame.color="black",cex=2,add=TRUE)
int_types=c('e','n','i','m')
color2D.matplot(xx,extremes=c(1:length(int_types)), border='white', axes=FALSE, xlab='', ylab='',main='',cellcolors=xx2)
legend(x=-15,y=num_play,legend=int_types,pch=22,pt.bg=pal,xpd=TRUE, bty='n')
dev.off()
"""















#Christmas tree orientation for food webs
R"""
library(ellipse)
maxTL = 8
TLcoords = list()
for (i in 1:maxTL) {
    whole_ell = ellipse(matrix(c(0,0.0,0.0,0),2,2),c(2/(i+3),1))
    low = which(whole_ell[,2] < 1 - i*0.5)
    TLcoords[[i]] = whole_ell[low,]
    TLcoords[[i]][,2] = TLcoords[[i]][,2] + i + 2 + rnorm(length(TLcoords[[i]][,2]),0,0.1)
    
}
plot(TLcoords[[1]],xlim=c(-2,2),ylim=c(0,maxTL+1))
for (i in 2:maxTL) {
    points(TLcoords[[i]])
}
"""

#0.4 if maxTL = 8
#0.2 if maxTL = 12