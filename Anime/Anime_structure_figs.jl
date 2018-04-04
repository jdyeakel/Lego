loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

d = load(string("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/structure.jld"));
#This loads the dictionary
rich = d["rich"];
sprich = d["sprich"];
turnover = d["turnover"];
mres_overlap = d["mres_overlap"];
conn = d["conn"];
prim_ext = d["prim_ext"];
sec_ext = d["sec_ext"];
status = d["status"];
lpot_col = d["lpot_col"];
degrees = d["degrees"];
trophic = d["trophic"];
avgdegree = d["avgdegree"];

reps = size(rich)[1];
tmax = size(rich)[2];

#trim zeros from degrees and trophic vectors
degrees_trim = Array{Array{Int64}}(reps);
trophic_trim = Array{Array{Float64}}(reps);
for i=1:reps
    degrees_trim[i] = degrees[i,find(!iszero,degrees[i,:])];
    trophic_trim[i] = trophic[i,find(!iszero,trophic[i,:])];
end

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
windowsize = 10;
#Connectance
R"""
library(RColorBrewer)
par(mfrow=c(3,1))
plot($(movingaverage(conn[1,2:tmax],windowsize)),log='x',type='l',ylim=c(0.02,0.06))
"""
for i=2:reps
    R"""lines($(movingaverage(conn[i,2:tmax],windowsize)))"""
end
#Specialization
R"""plot($(movingaverage(avgdegree[1,2:tmax],windowsize)),log='x',type='l',ylim=c(0,8))"""
for i=2:reps
    R"""lines($(movingaverage(avgdegree[i,2:tmax],windowsize)))"""
end
#Trophic Overlap
R"""plot($(movingaverage(mres_overlap[1,2:tmax],windowsize)),log='x',type='l',ylim=c(0,0.1))"""
for i=2:reps
    R"""lines($(movingaverage(mres_overlap[i,2:tmax],windowsize)))"""
end


m_conn = mapslices(median,conn,1);
m_mres_overlap = mapslices(median,mres_overlap,1);
R"""
library(RColorBrewer)
par(mfrow=c(2,1))
plot(seq(1,$tmax),$(m_conn),log='x',pch='.',ylim=c(0,0.05))
plot(seq(1,$tmax),$(m_mres_overlap),log='x',pch='.',ylim=c(0,0.05))
"""
