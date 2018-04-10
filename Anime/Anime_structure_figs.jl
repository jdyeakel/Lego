loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

d = load(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"));
#This loads the dictionary
rich = d["rich"];
sprich = d["sprich"];
turnover = d["turnover"];
mres_overlap = d["mres_overlap"];
conn = d["conn"];
conn_ind = d["conn_ind"];
prim_ext = d["prim_ext"];
sec_ext = d["sec_ext"];
status = d["status"];
lpot_col = d["lpot_col"];
avgdegree = d["avgdegree"];
degrees = d["degrees"];
trophic = d["trophic"];

reps = size(rich)[1];
tmax = size(rich)[2];
S = size(degrees)[3];



################
#Connectance
################

#There is a trail of zeros up front.
#Trim the zeros and fit curve
R"M_c=list(); tic=0";
parms_c = Array{Float64}(0,3);
for i=1:reps
    conn_trim = conn[i,find(!iszero,conn[i,:])];
    x = collect(1:length(conn_trim));
    y = conn_trim;
    R"""
    y <- $y
    x <- $x
    cf <- c(0,0,0)
    m <- try(nls(y ~ a + b * I(x^z), start = list(a = 1, b = 1, z = -1)),silent=TRUE)
        if (class(m) != "try-error") {
            tic = tic + 1;
            M_c[[tic]]=m;
            cf = coef(M_c[[tic]])
            }
    """
    # push!(b,Float64(R"coef(M[[tic]])[2]"));
    parms_c = vcat(parms_c,@rget(cf)')
end

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/conn_time.pdf");
R"""
pdf($namespace,height=5,width=6)
plot(fitted(M_c[[1]]),type='l', col = '#4292E515', lwd = 2,log='x',ylim=c(0,0.25),
xlab='Time',ylab='Connectance')
"""
for i=2:length(@rget(M_c))
    R"""
    lines(fitted(M_c[[$i]]), col = '#4292E515', lwd = 2)
    """
end
R"dev.off()"

#Percent of connectance trajectories that start high and asymptote to lower values
length(find(x->x>0,parms_c[:,2]))/reps

i=104
R"M_c[[$i]]"
R"""
plot(fitted(M_c[[$i]]),type='l', col = '#80808050', lwd = 2,log='x',ylim=c(0,0.2))
points($(conn[i,find(!iszero,conn[i,:])]))
"""

conn_trim = Array{Float64}(reps,tmax)*0;
conn_ind_trim = Array{Float64}(reps,tmax)*0;
for i=1:reps
    conn_rm = conn[i,find(!iszero,conn[i,:])];
    conn_ind_rm = conn_ind[i,find(!iszero,conn[i,:])];
    conn_trim[i,1:length(conn_rm)] = conn_rm;
    conn_ind_trim[i,1:length(conn_ind_rm)] = conn_ind_rm;
end
seq = [2;5;10;50;100;1000;2000];

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/conn_time2.pdf");
R"""
pdf($namespace,height=5,width=6)
boxplot($(conn_trim[:,seq]),ylim=c(0,0.1),outline=FALSE,names=$seq,
xlab='Time',ylab='Connectance',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(vec(mapslices(mean,conn_trim[:,seq],1))),ylim=c(0,0.1),pch=16)
lines($(vec(mapslices(mean,conn_trim[:,seq],1))),ylim=c(0,0.1),lwd=2)
dev.off()
"""

#Indirect connectance
namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/connind_time2.pdf");
R"""
pdf($namespace,height=5,width=6)
boxplot($(conn_ind_trim[:,seq]),ylim=c(0,0.1),outline=FALSE,names=$seq,
xlab='Time',ylab='Connectance',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(vec(mapslices(mean,conn_ind_trim[:,seq],1))),ylim=c(0,0.1),pch=16)
lines($(vec(mapslices(mean,conn_ind_trim[:,seq],1))),ylim=c(0,0.1),lwd=2)
dev.off()
"""



##############################
#Trophic Overlap
##############################

#Trophic Overlap
#Seems to be best summarized with boxplot
mres_overlap_trim = Array{Float64}(reps,tmax)*0;
for i=1:reps
    mres_overlap_rm = mres_overlap[i,find(x->x>0,mres_overlap[i,:])];
    mres_overlap_trim[i,1:length(mres_overlap_rm)] = mres_overlap_rm;
end
seq = [2;5;10;50;100;1000;2000];

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/trophicoverlap_time.pdf");
R"""
pdf($namespace,height=5,width=6)
boxplot($(mres_overlap_trim[:,seq]),ylim=c(0,0.1),outline=FALSE,names=$seq,
xlab='Time',ylab='Trophic overlap',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(vec(mapslices(mean,mres_overlap_trim[:,seq],1))),ylim=c(0,0.1),pch=16)
lines($(vec(mapslices(mean,mres_overlap_trim[:,seq],1))),ylim=c(0,0.1),lwd=2)
dev.off()
"""

mres_overlap_trim = mres_overlap[1,find(x->x>0,mres_overlap[1,:])];
R"""plot($(mres_overlap_trim),log='x',type='l',ylim=c(0,0.1))"""
for i=2:reps
    mres_overlap_trim = mres_overlap[i,find(x->x>0,mres_overlap[i,:])];
    R"""lines($(mres_overlap_trim))"""
end





##############################
#Specialization/Generalization
##############################
#Calculate from degrees
dt = Array{Float64}(reps,tmax);
for i=1:reps
    dvec = degrees[i,:,:];
    for t = 1:tmax
        dt[i,t] = mean(dvec[t,find(!iszero,dvec[t,:])]);
    end
end
meandegree = mapslices(mean,dt,1);

seq = [2;5;10;50;100;1000;2000];
namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/avgdegree_time.pdf");
R"""
pdf($namespace,height=5,width=6)
boxplot($(dt[:,seq]),ylim=c(0,30),outline=FALSE,names=$seq,
xlab='Time',ylab='Realized mean degree',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(vec(mapslices(mean,dt[:,seq],1))),ylim=c(0,0.1),pch=16)
lines($(vec(mapslices(mean,dt[:,seq],1))),ylim=c(0,0.1),lwd=2)
dev.off()
"""





# Somehow display how the *realized* degree distribution changes
seq = [10;50;100;1000;2000];

mdegt = Array{Float64}(length(seq),S)*0;
t_ic = 0;
for t = seq
    t_ic = t_ic + 1;
    deg = degrees[:,t,:];
    #Sort each row
    degsort = sort(deg,2,rev=true);
    #which column is the last non-zero?
    lastcol = find(iszero,sum(degsort,1))[1]-1;
    mdeg = Array{Float64}(lastcol);
    #Take means but ignore zeros for each column through lascol
    for i=1:lastcol
        mdeg[i] = median(degsort[!iszero.(degsort[:,i]),i]);
    end
    mdegt[t_ic,1:length(mdeg)]=mdeg;
end

namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/figures/degreedist_time.pdf");
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
pal = brewer.pal($(length(seq)),'Set1')
plot($(mdegt[1,!iszero.(mdegt[1,:])]),xlim=c(1,200),ylim=c(1,50),log='y',col=pal[1],type='l',lwd=2,xlab = 'Number of species', ylab='Median degree')
"""
for i=2:length(seq)
    R"""
    lines($(mdegt[i,!iszero.(mdegt[i,:])]),col=pal[$i],lwd=2)
    """
end
R"dev.off()"





i=1;
t=800;
deg = reverse(sort(degrees[i,t,find(!iszero,degrees[i,t,:])]));
R"""
plot($deg,ylim=c(1,50),log='y')
"""














###############
# OLD
###############


#Connectance and trophic width
R"""
library(RColorBrewer)
par(mfrow=c(2,1))
plot(seq(2,$tmax),$(conn[1,2:tmax]),log='x',type='l',ylim=c(0.0,0.5))
"""
for i=2:reps
    R"""lines(seq(2,$tmax),$(conn[i,2:tmax]))"""
end



R"""plot(seq(2,$tmax),$(mres_overlap[1,2:tmax]),log='x',type='l',ylim=c(0.0,0.25))"""
for i=2:reps
    R"""lines(seq(2,$tmax),$(mres_overlap[i,2:tmax]))"""
end

#Connectance, specialization, Trophic Overlap ~ moving average
windowsize = 1;
conn_ma = movingaverage(conn[:,2:tmax]',windowsize)';
avgdegree_ma = movingaverage(avgdegree[:,2:tmax]',windowsize)';
mres_overlap_ma = movingaverage(mres_overlap[:,2:tmax]',windowsize)';
tmax_ma = size(conn_ma)[2];
#Connectance
R"""
library(RColorBrewer)
par(mfrow=c(3,1))
plot($(conn_ma[1,2:tmax_ma]),log='x',type='l',ylim=c(0,0.1))
"""
for i=2:reps
    R"""lines($(conn_ma[i,2:tmax_ma]))"""
end
#Specialization
R"""plot($(avgdegree_ma[1,2:tmax_ma]),log='x',type='l',ylim=c(0,30))"""
for i=2:reps
    R"""lines($(avgdegree_ma[i,2:tmax_ma]))"""
end
#Trophic Overlap
R"""plot($(mres_overlap_ma[1,2:tmax_ma]),log='x',type='l',ylim=c(0,0.1))"""
for i=2:reps
    R"""lines($(mres_overlap_ma[i,2:tmax_ma]))"""
end


m_conn = mapslices(median,conn,1);
m_mres_overlap = mapslices(median,mres_overlap,1);
R"""
library(RColorBrewer)
par(mfrow=c(2,1))
plot(seq(1,$tmax),$(m_conn),log='x',pch='.',ylim=c(0,0.05))
plot(seq(1,$tmax),$(m_mres_overlap),log='x',pch='.',ylim=c(0,0.05))
"""
