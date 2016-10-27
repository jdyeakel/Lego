library(igraph)
library(plotrix)
library(RColorBrewer)
library(animation)

t=1000;
comgen <- read.table(paste(path.expand("~"),"/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen_t",t,".csv",sep=""),header=FALSE);
image(as.matrix(comgen),col=rev(gray.colors(2)))

image(matrix(as.matrix(comgen),500,100),col=c("white","black"))


tmax = 1000;
ani.options(interval=.025)
saveGIF({
  for (t in 1:tmax) {
    rl <- as.matrix(read.table(paste(path.expand("~"),"/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen_t",t,".csv",sep=""),header=FALSE));
    image(rl,col=rev(gray.colors(2)))
    #points(rwloc_frame[[i]],col="green")
  }
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_Lego/Lego_Program/movies/Lego.gif");
