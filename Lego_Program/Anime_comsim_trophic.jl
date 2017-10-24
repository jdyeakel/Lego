
@everywhere using Distributions
#@everywhere using Gadfly
@everywhere using RCall
@everywhere using HDF5
@everywhere using JLD


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsim.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/repsimintsingle.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_species.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/build_template_degrees.jl")


@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/initiate_comm_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/colonizesingle_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/extinct_func2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/sim_func.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/trophicalc2.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/trophicwidth.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/Jmatrix.jl")
@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/PSWebs.jl")

@everywhere include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/movingaverage.jl")

#Read-only variables
#Establish colonization and extinction rates
rate_col = 1;
#Establish thresholds
a_thresh = 0;

tmax = 500;
reps=50;


#Establish community template
S = 400;
probs = [
p_n=0.004,
p_a=0.01,
p_m=0.002,
p_i= 1 - sum([p_n,p_m,p_a]) #Ignore with 1 - pr(sum(other))
]


n_thresh = 0.2;
trophicload = 2;
ppweight = 1/4;

runpsw=true;

sprich,
rich,
conn,
#comgen,
ext_prim,
ext_sec,
tw,
twind,
psw,
pswind = repsimintsingle(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight,runpsw);


####################
# Jacobian stability
####################

R"""
library(RColorBrewer)
pal = brewer.pal(3,"Set1")
par(mfrow=c(2,1))
plot($(movingaverage(psw[:,1],10)),col=pal[1],type='l',ylim=c(0,1))
lines($(movingaverage(pswind[:,1],10)),col=pal[2])
"""
for i=2:reps
    R"""
    lines($(movingaverage(psw[:,i],10)),col=pal[1],type='l',ylim=c(0,1))
    lines($(movingaverage(pswind[:,i],10)),col=pal[2])
    """
end
R"""
plot($(sprich[:,1]),$(psw[:,1]),col=pal[1],pch=16,cex=0.8)
points($(sprich[:,1]),$(pswind[:,1]),col=pal[2],pch=16,cex=0.8)
"""
for i=2:reps
    R"""
    points($(sprich[:,i]),$(psw[:,i]),col=pal[1],pch=16,cex=0.8)
    points($(sprich[:,i]),$(pswind[:,i]),col=pal[2],pch=16,cex=0.8)
    """
end


#PSW at time t vs. Primary extinctions at t+1
delay=1;
R"""
plot($(pswind[:,1][1:tmax-delay]),$(ext_prim[:,1][(delay+1):tmax]),ylab='# primary extinctions @ time t+1',xlab='PSW @ time t',col=pal[1],pch=16,cex=0.4,ylim=c(0,3))
"""
for i=2:reps
    R"""
    points($(pswind[:,i][1:tmax-delay]),jitter($(ext_prim[:,i][(delay+1):tmax])),ylab='Number prim extinctions @ time t+1',xlab='PSW @ time t',col=pal[1],pch=16,cex=0.4)
    """
end
for i=0:3
    keep = find(x->x==i,ext_prim);
    psw_ext = pswind[keep];
    println(median(psw_ext))
    if i==0
        R"""
        library(RColorBrewer)
        pal = brewer.pal(3,"Set1")
        hist($psw_ext,col=pal[$i+1],freq=FALSE,breaks=20)
        """
    else
        R"""hist($psw_ext,add=TRUE,col=pal[$i+1],freq=FALSE,breaks=20)"""
    end
end



i=4
R"plot($(tw[:,i]),type='l')"


R"plot(apply($(tw),1,median),type='l',log='x',xlab='Time',ylab='Trophic overlap mean')"



timeatmax = zeros(Int64,tmax);
timeatmax_ind = zeros(Int64,tmax);
tw_aligned = zeros(tmax+1,reps);
twind_aligned = zeros(tmax+1,reps);
time = collect(1:tmax);
firstjump = find(x->x>0,tw[:,1])[1];
firstjump_ind = find(x->x>0,tw[:,1])[1];
timeatmax[1] = find(x->x==maximum(tw[:,1]),tw[:,1])[1];
timeatmax_ind[1] = find(x->x==maximum(twind[:,1]),twind[:,1])[1];
aseq = collect(firstjump:tmax);
aseq_ind = collect(firstjump_ind:tmax);
tw_aligned[2:length(aseq)+1,1] = tw[aseq,1];
twind_aligned[2:length(aseq_ind)+1,1] = twind[aseq_ind,1];
for i=2:reps
    firstjump = find(x->x>0,tw[:,i])[1];
    firstjump_ind = find(x->x>0,twind[:,i])[1];
    timeatmax[i] = find(x->x==maximum(tw[:,i]),tw[:,i])[1];
    timeatmax_ind[i] = find(x->x==maximum(twind[:,i]),twind[:,i])[1];
    aseq = collect(firstjump:tmax);
    aseq_ind = collect(firstjump_ind:tmax);
    tw_aligned[2:length(aseq)+1,i] = tw[aseq,i];
    twind_aligned[2:length(aseq_ind)+1,i] = twind[aseq_ind,i];
end

R"""
par(mfrow=c(1,1));
plot($(time-firstjump+1),$(tw[:,1]),type='l',xlab='Time',ylab='Trophic overlap mean',ylim=c(0,0.2),xlim=c(0,$tmax))
"""
for i=2:reps
  R"""
  lines($(time-firstjump+1),$(tw[:,i]))
  """
end

namespace = "$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/figures/fig_trophicwidth.pdf"
R"""
library(RColorBrewer)
library(boot)
pal = brewer.pal(3,'Set1')
mtw = apply($tw_aligned,1,mean);
mtwind = apply($twind_aligned,1,mean);
sdtw = apply($tw_aligned,1,sd);
plot(mtw,xlim=c(1,$tmax-50),type='l',col=pal[1],lwd=2,ylab='Trophic overlap',xlab='Time after growth',ylim=c(0,max(cbind(mtw,mtwind))))
lines(mtwind,xlim=c(1,$tmax-50),type='l',col=pal[1],lwd=2,lty=2)
#lines(mtw-sdtw,lty=2)
#lines(mtw+sdtw,lty=2)
#polygon(x=c(seq(1,length(mtw)),rev(seq(1,length(mtw)))),y=c(mtw+sdtw,rev(mtw-sdtw)),col=paste(pal[1],50,sep=''))
#lines(mtw,xlim=c(1,$tmax-50),type='l',col=pal[1],lwd=2)
"""

#Plot the time at which trophic overlap is maximized
R"""
hist($timeatmax_ind,breaks=20,col='gray',xlab='Time at trophic overlap maximum',ylim=c(0,80))
hist($timeatmax,breaks=20,col='lightblue',add=TRUE)
"""



