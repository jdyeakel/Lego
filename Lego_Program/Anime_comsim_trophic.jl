include("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/src/loadfuncs.jl")

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

save("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_psw/psw.jld",
"sprich",sprich,
"rich",rich,
"conn",conn,
#comgen,
"ext_prim",ext_prim,
"ext_sec",ext_sec,
"tw",tw,
"twind",twind,
"psw",psw,
"pswind",pswind);

d=load("$(homedir())/Dropbox/PostDoc/2014_Lego/Lego_Program/data/comsim_psw/psw.jld")
sprich=d["sprich"];
rich=d["rich"];
conn=d["conn"];
ext_prim=d["ext_prim"];
ext_sec=d["ext_sec"];
tw = d["tw"];
twind=d["twind"];
psw=d["psw"];
pswind=d["pswind"];


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

#PSW vs. extinctions
R"""
op=10;
plot($(sprich[:,1]),$(psw[:,1]),col=paste(pal[1],op,sep=''),pch=16,cex=0.8,ylim=c(0,1),xlim=c(0,max($sprich)))
points($(sprich[:,1]),$(pswind[:,1]),col=paste(pal[2],op,sep=''),pch=16,cex=0.8)
"""
for i=2:reps
    R"""
    points($(sprich[:,i]),$(psw[:,i]),col=paste(pal[1],op,sep=''),pch=16,cex=0.8)
    points($(sprich[:,i]),$(pswind[:,i]),col=paste(pal[2],op,sep=''),pch=16,cex=0.8)
    """
end


#PSW at time t vs. Primary extinctions at t+1
delay=1;
R"""
plot($(pswind[:,1][1:tmax-delay]),$(ext_prim[:,1][(delay+1):tmax]),ylab='# primary extinctions @ time t+1',xlab='PSW @ time t',col=pal[2],pch=16,cex=0.4,ylim=c(0,5))
"""
for i=2:reps
    R"""
    points($(pswind[:,i][1:tmax-delay]),jitter($(ext_prim[:,i][(delay+1):tmax])),col=pal[2],pch=16,cex=0.4)
    """
end


R"psw_extlist = list()"
R"sprich_extlist = list()"
# R"diffpsw_extlist = list()"
delay=1;
ext_prim_alt = ext_prim[delay+1:tmax,:];
sprich_alt = sprich[1:tmax-delay,:];
psw_alt = psw[1:tmax-delay,:];
pswind_alt = pswind[1:tmax-delay,:];
for i=0:5
    keep = find(x->x==i,ext_prim_alt);
    psw_ext = pswind_alt[keep];
    # diffpsw = psw[2:tmax,:] - psw[1:tmax-1,:];
    # diffpsw_ext = diffpsw[keep];
    sprich_ext = sprich_alt[keep];
    println(median(psw_ext))
    R"psw_extlist[[$i+1]] = $psw_ext;"
    R"sprich_extlist[[$i+1]] = $sprich_ext;"
    # if i==0
    #     R"""
    #     library(RColorBrewer)
    #     pal = brewer.pal(3,"Set1")
    #     hist($psw_ext,col=pal[$i+1],freq=FALSE,breaks=20)
    #     """
    # else
    #     R"""hist($psw_ext,add=TRUE,col=pal[$i+1],freq=FALSE,breaks=20)"""
    # end
end
R"boxplot(psw_extlist,names=c('0','1','2','3','4','5'),ylab='PSW @ time t',xlab='# of extinctions @ time t+1',ylim=c(0,1))";
R"boxplot(sprich_extlist,names=c('0','1','2','3','4','5'),ylab='SpRich @ time t',xlab='# of extinctions @ time t+1')";

r=1
diffpsw = psw[2:tmax,r] - psw[1:tmax-1,r];
R"plot($(ext_prim[1:tmax-1,r]),$diffpsw)"
for r = 1:reps
    diffpsw = psw[2:tmax,r] - psw[1:tmax-1,r];
    R"points($(ext_prim[1:tmax-1,r]),$diffpsw)"
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



