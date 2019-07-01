

if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end


#Search for parameters that match niche model
# annealtime = 500;
# 
# species = Array{Float64}(undef,annealtime,2);
# conn = Array{Float64}(undef,annealtime,2);
# mdegree = Array{Float64}(undef,annealtime,2);
# stdindegree = Array{Float64}(undef,annealtime,2);
# stdoutdegree = Array{Float64}(undef,annealtime,2);
#starting error
# global errvec_old = [5];
#starting temperature


cn = 3.14;
ce = 1.41;
cp = 0; #1.0;
p_n = 0.002;
p_a = 0.01;

reps = 1000;
# 
# enig_species = Array{Float64}(undef,reps);
# enig_conn = Array{Float64}(undef,reps);
# enig_mdegree = Array{Float64}(undef,reps);
# enig_stdindegree = Array{Float64}(undef,reps);
# enig_stdoutdegree = Array{Float64}(undef,reps);
# 
# niche_species = Array{Float64}(undef,reps);
# niche_conn = Array{Float64}(undef,reps);
# niche_mdegree = Array{Float64}(undef,reps);
# niche_stdindegree = Array{Float64}(undef,reps);
# niche_stdoutdegree = Array{Float64}(undef,reps);

err_species = Array{Float64}(undef,reps);
err_conn = Array{Float64}(undef,reps);
err_mdegree = Array{Float64}(undef,reps);
err_stdindegree = Array{Float64}(undef,reps);
err_stdoutdegree = Array{Float64}(undef,reps);




S = 200;
maxits = 2000;
SOprobs = (
p_n=p_n,
p_a=p_a
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


#expected objects per species
lambda = 0;
athresh = 0;
nthresh = 1.0;
MaxN = convert(Int64,floor(S + S*lambda));


enig_species = SharedArray{Float64}(reps);
enig_conn = SharedArray{Float64}(reps);
enig_mdegree = SharedArray{Float64}(reps);
enig_stdindegree = SharedArray{Float64}(reps);
enig_stdoutdegree = SharedArray{Float64}(reps);
@sync @distributed for j = 1:reps
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    sprich,rich,clock,CID,events = assembly(
        int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
        athresh,nthresh,maxits,cn,ce,cp);

    Aenigma = a_b[CID[:,maxits],CID[:,maxits]];
    Senigma = size(Aenigma)[1];
    Aenigma = Aenigma[2:Senigma,2:Senigma];
    Senigma -= 1;
    Cenigma = sum(Aenigma)/(Senigma^2);

    enig_species[j] = Senigma;
    enig_conn[j] = Cenigma;
    enig_mdegree[j] = mean(sum(Aenigma,dims=2));
    enig_stdindegree[j] = std(sum(Aenigma,dims=2));
    enig_stdoutdegree[j] = std(sum(Aenigma,dims=1));
end


# 
# enig_species[r] = mean(jspecies[vec(findall(!isnan,jspecies))]);
# enig_conn[r] = mean(jconn[vec(findall(!isnan,jconn))]);
# enig_mdegree[r] = mean(jmdegree[vec(findall(!isnan,jmdegree))]);
# enig_stdindegree[r] = mean(jstdindegree[vec(findall(!isnan,jstdindegree))]);
# enig_stdoutdegree[r] = mean(jstdoutdegree[vec(findall(!isnan,jstdoutdegree))]);


#Simulate a bunch of nichemodelwebs
niche_species = SharedArray{Float64}(reps);
niche_conn = SharedArray{Float64}(reps);
niche_mdegree = SharedArray{Float64}(reps);
niche_stdindegree = SharedArray{Float64}(reps);
niche_stdoutdegree = SharedArray{Float64}(reps);

@sync @distributed for i=1:reps
    #NICHE MODEL COMPARISON
    #NOTE: C CANNOT BE EQUAL TO OR GREATER THAN 0.5
    #NOTE Violates beta parameters!
    Aniche, n = nichemodelweb(Int64(floor(enig_species[i]))+20,enig_conn[i]);

    #make measurements
    Aniche = Array{Int64}(Aniche);

    Sniche = size(Aniche)[1];
    Cniche = sum(Aniche)/Sniche^2;

    niche_species[i] = Sniche;
    #Connectance
    niche_conn[i] = Cniche;

    #mean degree
    niche_mdegree[i] = mean(sum(Aniche,dims=2));

    #indegree
    niche_stdindegree[i] = std(sum(Aniche,dims=2));

    #outdegree
    niche_stdoutdegree[i] = std(sum(Aniche,dims=1));
end

# # 
# #ignore NAN
# 
# niche_species[r] = mean(ispecies);
# 
# #Connectance
# niche_conn[r] = mean(iconn);
# 
# #mean degree
# niche_mdegree[r] = mean(imdegree);
# 
# 
# #indegree
# niche_stdindegree[r] = mean(istdindegree);
# 
# 
# #outdegree
# niche_stdoutdegree[r] = mean(istdoutdegree);

for r=1:reps
    #Calculate error
    err_species[r] = sqrt((enig_species[r] - niche_species[r])^2); #/mean(species[r,1]); #/std(ispecies);
    err_conn[r] = sqrt((enig_conn[r] - niche_conn[r])^2); #/mean(conn[r,1]); #/std(iconn);
    err_mdegree[r] = sqrt((enig_mdegree[r] - niche_mdegree[r])^2); #/mean(mdegree[r,1]); #/std(imdegree);
    err_stdindegree[r] = sqrt((enig_stdindegree[r] - niche_stdindegree[r])^2); #/mean(stdindegree[r,1]); #/std(istdindegree);
    err_stdoutdegree[r] = sqrt((enig_stdoutdegree[r] - niche_stdoutdegree[r])^2); #/mean(stdoutdegree[r,1]); #/std(istdoutdegree);
end

namespace = "data/niche/errordist_cp0.jld";
filename = smartpath(namespace);
@save filename enig_species enig_conn enig_mdegree enig_stdindegree enig_stdoutdegree niche_species niche_conn niche_mdegree niche_stdindegree niche_stdoutdegree err_species err_conn err_mdegree err_stdindegree err_stdoutdegree;

namespace = "data/niche/errordist_ce0.jld";
filename = smartpath(namespace);
@load filename enig_species enig_conn enig_mdegree enig_stdindegree enig_stdoutdegree niche_species niche_conn niche_mdegree niche_stdindegree niche_stdoutdegree err_species err_conn err_mdegree err_stdindegree err_stdoutdegree;

#Plot normalized error

filename = "figures/niche/error.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,height=12,width=3)
par(mfrow=c(5,1))
hist($(err_species),main="",xlab="Error (Species)",col='darkgray')
hist($(err_conn),main="",xlab="Error (Connectance)",col='darkgray')
hist($(err_mdegree),main="",xlab="Error (Mean degree)",col='darkgray')
hist($(err_stdindegree),main="",xlab="Error (In-Degree SD)",col='darkgray')
hist($(err_stdoutdegree),main="",xlab="Error (Out-Degree SD)",col='darkgray')
dev.off()
"""


filename = "figures/niche/normerror.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,height=12,width=3)
par(mfrow=c(5,1))
hist($(err_species ./ std(enig_species)),main="",xlab="Normalized Error (Species)",col='darkgray')
hist($(err_conn ./ std(enig_conn)),main="",xlab="Normalized Error (Connectance)",col='darkgray')
hist($(err_mdegree ./ std(enig_mdegree)),main="",xlab="Normalized Error (Mean degree)",col='darkgray')
hist($(err_stdindegree ./ std(enig_stdindegree)),main="",xlab="Normalized Error (In-Degree SD)",col='darkgray')
hist($(err_stdoutdegree ./ std(enig_stdoutdegree)),main="",xlab="Normalized Error (Out-Degree SD)",col='darkgray')
dev.off()
"""

filename = "figures/niche/errorscatter_cp0.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,height = 10,width=10)
par(mfrow=c(2,3))
plot($niche_species,$enig_species,pch='.',ylim=c(50,200),xlim=c(50,200),xlab="Niche (species)",ylab="ENIgMa (species)")
lines(seq(-10,500),seq(-10,500),lty=2)
plot($niche_conn,$enig_conn,pch='.',ylim=c(0,0.025),xlim=c(0,0.025),xlab="Niche (connectance)",ylab="ENIgMa (connectance)")
lines(seq(-10,10),seq(-10,10),lty=2)
plot($niche_mdegree,$enig_mdegree,pch='.',ylim=c(0,3),xlim=c(0,3),xlab="Niche (mean degree)",ylab="ENIgMa (mean degree)")
lines(seq(-10,10),seq(-10,10),lty=2)
plot($niche_stdindegree,$enig_stdindegree,pch='.',ylim=c(0.5,4),xlim=c(0.5,4),xlab="Niche (in-degree SD)",ylab="ENIgMa (in-degree SD)")
lines(seq(-10,10),seq(-10,10),lty=2)
plot($niche_stdoutdegree,$enig_stdoutdegree,pch='.',ylim=c(0.5,2.5),xlim=c(0.5,2.5),xlab="Niche (out-degree SD)",ylab="ENIgMa (out-degree SD)")
lines(seq(-10,10),seq(-10,10),lty=2)
dev.off()
"""

