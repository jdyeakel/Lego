library(igraph)
library(plotrix)
library(RColorBrewer)
library(animation)

t=1;
comgen <- read.table(paste(path.expand("~"),"/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen_t",t,".csv",sep=""),header=FALSE);
image(as.matrix(comgen),col=rev(gray.colors(2)))




tmax = 1000;
ani.options(interval=.025)
saveGIF({
  for (t in 1:tmax) {
    rl <- as.matrix(read.table(paste(path.expand("~"),"/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen_t",t,".csv",sep=""),header=FALSE));
    image(matrix(rl,500,100),col=c("white","black"))
    #points(rwloc_frame[[i]],col="green")
  }
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/movies/Lego.gif")
