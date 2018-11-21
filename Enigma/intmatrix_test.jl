if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
# @everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/assembly2.jl")


S = 40;

maxits = 1000;

SOprobs = (
p_n=0.01,
p_a=0.03
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


#expected objects per species
lambda = 0.0;
MaxN = convert(Int64,floor(S + S*lambda));

int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

a_b,n_b,i_b,m_b,n_b0,sp_v,int_id = preamble_defs(int_m);

#Reorganize to clump objects
objects = deleteat!(findall(x->x=='i',diag(int_m)),1);
objectssort = objects[sortperm(vec(sum(m_b[objects,objects],dims=1)),rev=true)];
# objectssort2 = objectssort[sortperm(vec(sum(m_b[objectssort,objectssort].+a_b[objectssort,objectssort].+n_b[objectssort,objectssort],2)),rev=false)];
species = setdiff(collect(1:length(diag(int_m))),objects);
speciessort = species[sortperm(vec(sum(a_b[species,species].+n_b[species,species],dims=1)),rev=true)];
speciessort2 = speciessort[sortperm(vec(sum(a_b[speciessort,speciessort].+n_b[speciessort,speciessort],dims=2)),rev=false)];
int_msort = int_m[[speciessort2;objectssort],[speciessort2;objectssort]];


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
pdf($namespace,height=5,width=6)
par(mar=c(1,1,1,4))
int_types=c('e','n','i','m')
color2D.matplot(xx,extremes=c(1:length(int_types)), border='white', axes=FALSE, xlab='', ylab='',main='',cellcolors=xx2)
legend(x=num_play+1,y=num_play,legend=int_types,pch=22,pt.bg=pal,xpd=TRUE, bty='n')
#text(x=rep(-0.8,length(objects)),y=num_play-objects+0.5,labels='o', xpd=TRUE,cex=0.6)
#text(x=objects-0.5,y=rep(num_play+0.8,length(objects)),labels='o', xpd=TRUE,cex=0.6)
text(x=-0.8,y=num_play-0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
text(x=0.8,y=num_play+0.8,labels=expression(paste('1',degree)), xpd=TRUE,cex=0.6)
dev.off()
"""

