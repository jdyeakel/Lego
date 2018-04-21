# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");
loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");

# namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/sim.jld");
namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/sim.jld");
d1 = load(namespace);
lpot_col = d1["lpot_col"];
status = d1["status"];
prim_ext = d1["prim_ext"];
sec_ext = d1["sec_ext"];
cid_r = d1["cid_r"];

reps = size(cid_r)[1];
tmax = size(cid_r)[2];
S = convert(Int64,size(cid_r)[3]/2);

# #We need to 'grab' individual cid_r's within the parallel loop
# #Export them here and then import them separately under parallel rep loop
# @time for r=1:reps
#     namespace = string("$(homedir())/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
#     # namespace = string("/$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
#     CID = cid_r[r,:,:];
#     save(namespace,
#     "CID", CID);
# end

#define cid_r everywhere


# NOTE: just analyze the timesteps necessary for the sequence

# seq = [10;25;50;100;200;1000;2000;];

#We want to catch the full assembly early on, and then skip ahead
seq = [collect(2:50);100;200;1000;2000;];
tseqmax = length(seq);

rich = SharedArray{Int64}(reps,tseqmax);
sprich = SharedArray{Int64}(reps,tseqmax);
turnover = SharedArray{Float64}(reps,tseqmax);
res_overlap = SharedArray{Float64}(reps,tseqmax);
conn = SharedArray{Float64}(reps,tseqmax);
conn_ind = SharedArray{Float64}(reps,tseqmax);
avgdegree = SharedArray{Float64}(reps,tseqmax);
res_overlap_dist = SharedArray{Float64}(reps,tseqmax,S);
degrees = SharedArray{Int64}(reps,tseqmax,S);
trophic = SharedArray{Float64}(reps,tseqmax,S);

@time @sync @parallel for r=1:reps
    #Read in the interaction matrix
    # namespace = string("$(homedir())/Dropbox/Postdoc/2014_Lego/Anime/data/simbasic/int_m",r,".jld");

    namespace_rep = string("$(homedir())/2014_Lego/Anime/data/simbasic/int_m",r,".jld");
    
    d2 = load(namespace_rep);
    int_m = d2["int_m"];
    sp_m = d2["sp_m"];
    t_m = d2["t_m"];
    tp_m = d2["tp_m"];
    tind_m = d2["tind_m"];
    mp_m = d2["mp_m"];
    mind_m = d2["mind_m"];
    
    namespace_cid = string("$(homedir())/2014_Lego/Anime/data/simbasic/cid_",r,".jld");
    d3 = load(namespace_cid);
    CID = d3["CID"];
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    #Analysis
    for t = 1:tseqmax
        
        #construct
        tstep = seq[t];
        cid = find(isodd,CID[tstep,:]);
        cid_old = find(isodd,CID[tstep-1,:]); #because we have this, seq can't start at t=1;
        
        rich[r,t], sprich[r,t], turnover[r,t], res_overlap[r,t], res_overlap_all, conn[r,t], conn_ind[r,t] = dynstructure(cid,cid_old,sp_v,a_b,tp_m,tind_m);     
        
        res_overlap_dist[r,t,1:length(res_overlap_all)] = res_overlap_all; 
        
        deg,troph = structure(S,cid,sp_v,tind_m);
        
        degrees[r,t,1:length(deg)] = deg;
        trophic[r,t,1:length(troph)] = troph;
        avgdegree[r,t] = mean(degrees[r,t,1:length(deg)]);
        
        
    end

end

h = load(string("$(homedir())/2014_Lego/Anime/data/intm_structure.jld"));
Pconn = h["Pconn"];
Pconn_ind = h["Pconn_ind"];
Pres_overlap_dist = h["Pres_overlap_dist"];
Pdegrees = h["Pdegrees"];
Ptl = h["Ptl"];


#Connectance
lfseq = find(x->x>1,diff(seq))[1];
#This is the initial assembly process
init_conn = conn[:,1:lfseq];
init_conn_ind = conn_ind[:,1:lfseq];
init_conn_trim = Array{Float64}(reps,lfseq)*0;
init_conn_ind_trim = Array{Float64}(reps,lfseq)*0;
for r=1:reps
    init_conn_rm = init_conn[r,find(!iszero,init_conn[r,:])];
    init_conn_ind_rm = init_conn_ind[r,find(!iszero,init_conn_ind[r,:])];
    init_conn_trim[r,1:length(init_conn_rm)] = init_conn_rm;
    init_conn_ind_trim[r,1:length(init_conn_ind_rm)] = init_conn_ind_rm;
end
bins = [2;5;10;50;100;1000;2000];
initsteps = bins[bins.<lfseq]; #use these locations for init
laststeps = bins[bins.>=lfseq]; #use these locations for the rest
lastbins = indexin(laststeps,seq);

#Stitch together
seq_stitch = [initsteps;laststeps];
conn_stitch = [init_conn_trim[:,initsteps] conn[:,lastbins]];

namespace = string("$(homedir())/conn_time2NEW.pdf");
R"""
pdf($namespace,height=5,width=6)
boxplot($(conn_stitch),ylim=c(0,0.1),outline=FALSE,names=$(seq_stitch),
xlab='Time',ylab='Connectance',
pars = list(boxwex = 0.4, staplewex = 0.5, outwex = 0.5),col='lightgray')
points($(vec(mapslices(mean,conn_stitch,1))),ylim=c(0,0.1),pch=16)
lines($(vec(mapslices(mean,conn_stitch,1))),ylim=c(0,0.1),lwd=2)
lines(seq(0.001,3000),rep(mean($Pconn),3000),lty=2)
dev.off()
"""



