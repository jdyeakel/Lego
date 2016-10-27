library(igraph)
library(plotrix)
library(RColorBrewer)


comgen <- read.table(paste(path.expand("~"),"/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comgen.csv",sep=""),header=FALSE);
image(as.matrix(comgen))


color2D.matplot(as.matrix(comgen), border="white", axes=FALSE, xlab="", ylab="",main="")



#Visualize matrix:
plot_matrix<-function(int_m, num_play){
xx=matrix(as.numeric(as.factor(int_m)),c(num_play,num_play))
par(mar=c(1,1,1,4))
int_types=levels(as.factor(int_m))
color2D.matplot(xx,extremes=c(1:length(int_types)), border="white", axes=FALSE, xlab="", ylab="",main="")
legend(x=num_play+1,y=num_play,legend=int_types,pch=22,pt.bg=c(1:5),xpd=TRUE, bty="n")
}
plot_matrix(comgen, )







ani.options(interval=.025)
saveGIF({
  for (i in 1:tmax) {
    rl <- r_frame[i,]
    image(matrix(rl,(L+2),(L+2)),col=c("white","black"))
    #points(rwloc_frame[[i]],col="green")
  }
},movie.name = "/Users/justinyeakel/Dropbox/PostDoc/2014_DiffusingForager/animations/resourceL50_pulse.gif")
