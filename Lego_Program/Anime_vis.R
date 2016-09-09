library(igraph)
library(RColorBrewer)
pal <- brewer.pal(3,"Set1")


fw <- read.table("/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/data_template/fw.csv",header=FALSE)

fw_g <- graph.adjacency(as.matrix(fw))
basal_pos <- 1
num_play = dim(fw)[1]
trophic <- sapply(1:vcount(fw_g),function(x){mean(shortest.paths(fw_g,basal_pos,which(fw[x,]==1)))+0})
#trophic[which(trophic==Inf)] <- 0
trophic[which(trophic=="NaN")] <- 0
coords <- cbind(runif(vcount(fw_g)),trophic); coords[basal_pos,] <- c(0.5,trophic[basal_pos])
par(mar=c(1,1,1,1))
plot(fw_g,layout=coords,vertex.size=5,edge.arrow.size=0.5,
     main=ecount(fw_g)/(num_play^2),vertex.label=NA,
     vertex.color=pal[2])
