if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end

filename = "data/steadystate/sim_settings.jld";
namespace = smartpath(filename);
# namespace = "$(homedir())/2014_Lego/Enigma/data/steadystate/sim_settings.jld";
@load namespace reps S maxits athresh nthresh lambda SSprobs SOprobs OOprobs;

# 
# d1 = load(namespace);
# reps = d1["reps"];
# S = d1["S"];
# maxits = d1["maxits"];
# athresh = d1["athresh"];
# nthresh = d1["nthresh"];


seq = [collect(2:50);100;200;500;1000;2000;4000];
tseqmax = length(seq);

rich = SharedArray{Int64}(reps,tseqmax);
sprich = SharedArray{Int64}(reps,tseqmax);
sprichinweb = SharedArray{Int64}(reps,tseqmax);
sprichinwebnoprim = SharedArray{Int64}(reps,tseqmax);

turnover = SharedArray{Float64}(reps,tseqmax);
res_overlap = SharedArray{Float64}(reps,tseqmax);
user_overlap = SharedArray{Float64}(reps,tseqmax);
conn = SharedArray{Float64}(reps,tseqmax);
conn_ind = SharedArray{Float64}(reps,tseqmax);
mutconn = SharedArray{Float64}(reps,tseqmax);
mutconn_ind = SharedArray{Float64}(reps,tseqmax);
avgdegree = SharedArray{Float64}(reps,tseqmax);
pc = SharedArray{Int64}(reps,tseqmax);
res_overlap_dist = SharedArray{Float64}(reps,tseqmax,S);
user_overlap_dist = SharedArray{Float64}(reps,tseqmax,S);
degrees = SharedArray{Int64}(reps,tseqmax,S);
trophic = SharedArray{Float64}(reps,tseqmax,S);

realizeddegree = SharedArray{Float64}(reps,tseqmax,S);
potentialdegree = SharedArray{Float64}(reps,tseqmax,S);

realizedG = SharedArray{Float64}(reps,tseqmax,S);
realizedGnoprim = SharedArray{Float64}(reps,tseqmax,S);
realizedGavgc = SharedArray{Float64}(reps,tseqmax,S);
potentialGavgc = SharedArray{Float64}(reps,tseqmax,S);

potentialG = SharedArray{Float64}(reps,tseqmax,S);
C = SharedArray{Float64}(reps,tseqmax);

#Save richness for colonizing phase
coltraj = SharedArray{Int64}(reps,1000);
clocktraj = SharedArray{Float64}(reps,1000);


#FIRST Calculate SSL
SS_linkdensity = SharedArray{Float64}(reps);
SS_indlinkdensity = SharedArray{Float64}(reps);
@sync @distributed for r=1:reps
    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;
    
    cid = findall(isodd,CID[:,4000]);
    spcid = intersect(sp_v,cid);
    
    cidinwebnoprim = sort(cid[findall(!iszero,vec(sum(a_b[cid,cid],dims=2)))]);
    
    #ONLY COUNT CONNECTED SPECIES AND NOT PURE PRIMARY PRODUCERS
    speciesatsteadystate = length(cidinwebnoprim);
    linksatsteadystate = sum(a_b[cidinwebnoprim,cidinwebnoprim]);
    indlinksatsteadystate = sum(tind_m[cidinwebnoprim,cidinwebnoprim]);
    
    # speciesatsteadystate = length(spcid);
    # linksatsteadystate = sum(a_b[spcid,spcid]);
    # indlinksatsteadystate = sum(tind_m[spcid,spcid]);
    
    SS_linkdensity[r] = linksatsteadystate/speciesatsteadystate;
    SS_indlinkdensity[r] = indlinksatsteadystate/speciesatsteadystate;
end

SSL = mean(SS_linkdensity);

@sync @distributed for r=1:reps
    #Read in the interaction matrix
    
    # if homedir() == "/home/z840"
    #     namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    # else
    #     namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/steadystate/int_m",r,".jld");
    # end
    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    
    # @load namespace_rep;
    # 
    # 
    # d2 = load(namespace_rep);
    # int_m = d2["int_m"];
    # tp_m = d2["tp_m"];
    # tind_m = d2["tind_m"];
    # mp_m = d2["mp_m"];
    # mind_m = d2["mind_m"];
    # 
    # if homedir() == "/home/z840"
    #     namespace = string("$(homedir())/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    # else
    #     namespace = string("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/data/steadystate/cid_",r,".jld");
    # end
    
    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;
    
    # d3 = load(namespace_cid);
    # CID = d3["CID"];
    # 
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
    
    #Save FULL colonizing phase
    coltraj[r,:] = sum(CID[:,1:1000],dims=1);
    clocktraj[r,:] = clock[1:1000];

    #Analysis
    for t = 1:tseqmax
        
        #construct
        tstep = seq[t];
        cid = findall(isodd,CID[:,tstep]);
        cid_old = findall(isodd,CID[:,tstep-1]); #because we have this, seq can't start at t=1;
        
        rich[r,t], sprich[r,t], turnover[r,t], res_overlap[r,t], user_overlap[r,t], res_overlap_all, user_overlap_all, conn[r,t], conn_ind[r,t], mutconn[r,t], mutconn_ind[r,t], pc_cid = dynstructure(cid,cid_old,sp_v,a_b,n_b0,tp_m,tind_m,mp_m,mind_m,int_id,athresh,nthresh);     
        
        pc[r,t] = length(pc_cid);
        
        res_overlap_dist[r,t,1:length(res_overlap_all)] = res_overlap_all;
        #Only save species user-overlap, and not object user-overlap
        user_overlap_dist[r,t,1:length(user_overlap_all)] = user_overlap_all; 
        
        deg,troph = structure(S,cid,sp_v,tind_m);
        
        degrees[r,t,1:length(deg)] = deg;
        trophic[r,t,1:length(troph)] = troph;
        avgdegree[r,t] = mean(degrees[r,t,1:length(deg)]);
        
        # realizeddegree[r,t] = mean(vec(sum(a_b[cid,[1;cid]],dims=2)));
        # potentialdegree[r,t] = mean(vec(sum(a_b[cid,:],dims=2)));
        realizeddegree[r,t,1:length(cid)] = vec(sum(a_b[cid,[1;cid]],dims=2));
        potentialdegree[r,t,1:length(cid)] = vec(sum(a_b[cid,:],dims=2));
        
        
        #NOTE species disconnected from web can be counted
        
        #these are species that are connected to the foodweb in any way
        cidinweb = sort(cid[findall(!iszero,vec(sum(a_b[cid,[1;cid]],dims=2)))]);
        
        sprichinweb[r,t] = length(cidinweb);
        
        C[r,t] = (sum(a_b[cidinweb,[1;cid]]))/((length(cidinweb))^2);
        
        realizedG[r,t,1:length(cidinweb)] = vec(sum(a_b[cidinweb,[1;cid]],dims=2)) .* (1/((sum(a_b[cidinweb,[1;cid]]))/(length(cidinweb))));
        
        potentialG[r,t,1:length(cidinweb)] = vec(sum(a_b[cidinweb,:],dims=2)) .* (1/((sum(a_b[cidinweb,[1;cid]]))/(length(cidinweb))));
        
        cidinwebnoprim = sort(cid[findall(!iszero,vec(sum(a_b[cid,cid],dims=2)))]);
        sprichinwebnoprim[r,t] = length(cidinwebnoprim);
        
        if sprichinwebnoprim[r,t] > 0
            realizedGnoprim[r,t,1:length(cidinwebnoprim)] = vec(sum(a_b[cidinwebnoprim,cid],dims=2)) .* (1/((sum(a_b[cidinwebnoprim,cid]))/(length(cidinwebnoprim))));
        end
        
        #compared to steady state C averaged across webs
        #AS CLOSE TO PIETCHNICK AS POSSIBLE
        # SSC = 0.0108;
        # SSL = 1.402144739810133;
        if sprichinwebnoprim[r,t] > 0
            # realizedGavgc[r,t,1:length(cidinweb)] = vec(sum(a_b[cidinweb,[1;cid]],dims=2)) .* (1/(SSC*length(cidinweb)));
            realizedGavgc[r,t,1:length(cidinwebnoprim)] = vec(sum(a_b[cidinwebnoprim,cid],dims=2)) .* (1/(SSL));
            potentialGavgc[r,t,1:length(cidinwebnoprim)] = vec(sum(a_b[cidinwebnoprim,2:S],dims=2)) .* (1/(SSL));
        end
        
    end

end


#Real degreea and trophic distribution

##################
#USE THESE 6/24/2019
##################

bins = [5;10;25;50;100;200;500;1000;2000;4000;];
# bins = [5;50;100;4000;];
seq2 = indexin(bins,seq);
tmaxdegree = zeros(Int64,length(seq2));
tmaxtrophic = zeros(Int64,length(seq2));
maxdegree = Array{Int64}(undef,reps,length(seq2));
maxtrophic = Array{Int64}(undef,reps,length(seq2));
freqdegreereps = Array{Array}(undef,reps);
freqtrophicreps = Array{Array}(undef,reps);
for r=1:reps
    freqdegreetime = Array{Array}(undef,length(seq2));
    freqtrophictime = Array{Array}(undef,length(seq2));
    let tic = 0
        for t=seq2
            tic += 1;
            alldegrees = degrees[r,t,:][degrees[r,t,:] .> 0];
            alltrophic = trophic[r,t,:][trophic[r,t,:] .>= 0.9];
            maxdegree[r,tic] = maximum(alldegrees);
            maxtrophic[r,tic] = round(Int64,maximum(alltrophic))+1;
            freqdegree = Array{Float64}(undef,maxdegree[r,tic]);
            freqtrophic = Array{Float64}(undef,maxtrophic[r,tic]+1);
            for i=1:maxdegree[r,tic]
                freqdegree[i] = length(findall(x->x==i,alldegrees))/length(alldegrees);
            end
            for i=0:maxtrophic[r,tic]
                freqtrophic[i+1] = length(findall(x->(x>=i && x<i+1),alltrophic))/length(alltrophic);
            end
            freqdegreetime[tic] = freqdegree;
            freqtrophictime[tic] = freqtrophic;
            tmaxdegree[tic] = maximum([maxdegree[r,tic],tmaxdegree[tic]])
            tmaxtrophic[tic] = maximum([maxtrophic[r,tic],tmaxtrophic[tic]])
        end
    end
    freqdegreereps[r] = freqdegreetime;
    freqtrophicreps[r] = freqtrophictime;
end
#lets ignore results with >20 trophic levels!
toignore = findall(!iszero,vec(sum(maxtrophic .> 20,dims=2)));
newreps = setdiff(collect(1:reps),toignore);


maxdegreeall = maximum(tmaxdegree);
# maxtrophic = maximum(tmaxtrophic)+1;
maxtrophicall = 21;
#Reorganizing
degreedistreps = zeros(Float64,length(newreps),length(seq2),maxdegreeall);
trophicdistreps = zeros(Float64,length(newreps),length(seq2),maxtrophicall);
let rtic=0
    for r=newreps
        rtic += 1;
        for t=1:length(seq2)
            degreedistreps[rtic,t,1:length(freqdegreereps[r][t])] = freqdegreereps[r][t];
            trophicdistreps[rtic,t,1:length(freqtrophicreps[r][t])] = freqtrophicreps[r][t];
        end
    end
end
#means over reps
meandegreedist = Array{Float64}(undef,length(seq2),maxdegreeall);
meantrophicdist = Array{Float64}(undef,length(seq2),maxtrophicall);
for t=1:length(seq2)
    meandegreedist[t,:] = mean(degreedistreps[:,t,:],dims=1);
    meantrophicdist[t,:] = mean(trophicdistreps[:,t,:],dims=1);
end


#Resource specialization of realized vs. potential niche over time

bins = [5;10;25;50;100;200;500;1000;2000;4000;];
# bins=seq;
seq2 = indexin(bins,seq);
meanpotentialdegree = Array{Float64}(undef,reps,length(seq2));
meanrealizeddegree = Array{Float64}(undef,reps,length(seq2));
rG = Array{Float64}(undef,reps,length(seq2));
pG = Array{Float64}(undef,reps,length(seq2));
propG = Array{Float64}(undef,reps,length(seq2));
propGavgc = Array{Float64}(undef,reps,length(seq2));
proppotGavgc = Array{Float64}(undef,reps,length(seq2));
propGnoprim = Array{Float64}(undef,reps,length(seq2));
uppertrophicpropG = Array{Float64}(undef,reps,length(seq2));
potpropG = Array{Float64}(undef,reps,length(seq2));
meantrophic = Array{Float64}(undef,reps,length(seq2));
meandegrees = Array{Float64}(undef,reps,length(seq2));

for r=1:reps
    for t=1:length(seq2)
        
        meanpotentialdegree[r,t] = mean(potentialdegree[r,seq2[t],:][potentialdegree[r,seq2[t],:].>0]);
        meanrealizeddegree[r,t] = mean(realizeddegree[r,seq2[t],:][realizeddegree[r,seq2[t],:].>0]);
        
        rG[r,t] = mean(realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0]);
        pG[r,t] = mean(potentialG[r,seq2[t],:][potentialG[r,seq2[t],:].>0]);
        
        # propG[r,t] = sum((realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0] .* C[r,seq2[t]]) .> (C[r,seq2[t]]))/sum(realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0]);
        
        # propG[r,t] = sum((realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0] ) .> (1))/sum(realizedG[r,seq2[t],:][realizedG[r,seq2[t],:].>0]);
        
        propG[r,t] = sum((realizedG[r,seq2[t],1:sprichinweb[r,seq2[t]]]) .> (1))/sprichinweb[r,seq2[t]];
        
        propGnoprim[r,t] = sum((realizedGnoprim[r,seq2[t],1:sprichinwebnoprim[r,seq2[t]]]) .> (1))/sprichinwebnoprim[r,seq2[t]];
        
        #USED IN MS
        #realized generality
        propGavgc[r,t] = sum((realizedGavgc[r,seq2[t],1:sprichinwebnoprim[r,seq2[t]]]) .> (1))/sprichinwebnoprim[r,seq2[t]];
        #potential generality
        proppotGavgc[r,t] = sum((potentialGavgc[r,seq2[t],1:sprichinwebnoprim[r,seq2[t]]]) .> (1))/sprichinwebnoprim[r,seq2[t]];
        
        potpropG[r,t] = sum((potentialG[r,seq2[t],1:sprichinweb[r,seq2[t]]]) .> (1))/sprichinweb[r,seq2[t]];
        
        uppertrophicpropG[r,t] = sum((realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])][realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])].>0]) .> (1))/sum(realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])][realizedG[r,seq2[t],findall(x->x<2,trophic[r,seq2[t],:])].>0]);
        
        meantrophic[r,t] = mean(trophic[r,seq2[t],:][trophic[r,seq2[t],:] .> 0]);
        meandegrees[r,t] = mean(degrees[r,seq2[t],:][degrees[r,seq2[t],:] .> 0]);
        
    end
end

#which reps have low values of propG durin timesteps 1:100?
earlygen = unique([findall(x->x<0.5,1 .- propG[:,3]);findall(x->x<0.5,1 .- propG[:,3])]);
earlysp = setdiff(collect(1:reps),earlygen);

mps = mean(meanpotentialdegree,dims=1);
mrs = mean(meanrealizeddegree,dims=1);


mpG = mean(potpropG,dims=1);
mrG = meanfinite(propG,1);
mrGgen = meanfinite(propG[earlygen,:],1);
mrGnoprim = meanfinite(propGnoprim,1);

mrGavgc = meanfinite(propGavgc,1);
mpGavgc = meanfinite(proppotGavgc,1);

minmrGnoprim = Array{Float64}(undef,length(seq2));
maxmrGnoprim = Array{Float64}(undef,length(seq2));
for i=1:length(seq2)
    minmrGnoprim[i]=quantile(propGnoprim[isnan.(propGnoprim[:,i]).==false,i],0.05);
    maxmrGnoprim[i]=quantile(propGnoprim[isnan.(propGnoprim[:,i]).==false,i],0.95);
end
# mrGuppertrophic = meanfinite(uppertrophicpropG,1);



# filename = "../manuscript/fig_trophic2.pdf";
# namespace = smartpath(filename);
filename = "figures/yog/fig_trophic2.pdf";
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pdf($namespace,height=5,width=6)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), 
   widths=c(1,1,1), heights=c(0.8,1,1))
par(oma = c(0.5, 1, 1, 1), mar = c(3, 4, 1, 1))
"""
rsample=sample(collect(1:reps),100);
r=rsample[1];
filename = "data/steadystate/int_m.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace int_m tp_m tind_m mp_m mind_m;
filename = "data/steadystate/cid.jld";
indices = [r];
namespace = smartpath(filename,indices);
@load namespace CID clock;
R"""
pal=brewer.pal(3,'Set1')
plot(seq(1,4000),$(vec(sum(CID,dims=1))),type='l',col=paste(pal[2],40,sep=''),xlim=c(0,4000),ylim=c(0,180),xlab='',ylab='',axes=FALSE)
axis(1)
axis(2,las=1)
title(ylab='Species richness', line=3, cex.lab=1.2)
title(xlab='Assembly time', line=2.0, cex.lab=1.2)
mtext(paste0("a"), side = 3, adj = -0.16, 
    line = 0.3,cex=1.2,font=2)
"""

for i=2:length(rsample)
    r=rsample[i];
    filename = "data/steadystate/int_m.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace int_m tp_m tind_m mp_m mind_m;
    filename = "data/steadystate/cid.jld";
    indices = [r];
    namespace = smartpath(filename,indices);
    @load namespace CID clock;
    # filename = "figures/yog/assembly_time.pdf"
    # namespace = smartpath(filename);
    R"""
    lines(seq(1,4000),$(vec(sum(CID,dims=1))),type='l',col=paste(pal[2],40,sep=''))
    """
end



#SPECIALIZATION
R"""
pal = brewer.pal($(length(seq2)),'Spectral')
timelabels = parse(text=c("5","10","25","50",paste("10","^2"),paste("2.10","^2"),paste("5*10","^2"),paste("10","^3"),paste("2*10","^3"),paste("4*10","^3")))
# par(mfrow=c(2,1))
plot($(seq[seq2]),1 - $mrGavgc, xlim=c(5,4000), ylim=c(0,1), xlab='', ylab='', log='x', cex.axis=0.85, pch=23, col='black', bg=pal, cex=1.5,axes=FALSE)
axis(2,at=seq(0,1,by=0.2),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=$(seq[seq2]),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Proportion specialists', line=3, cex.lab=1.2)
title(xlab='Assembly time', line=2.0, cex.lab=1.2)
mtext(paste0("b"), side = 3, adj = -0.4, 
    line = -0.6,cex=1.2,font=2)
"""
for i=1:length(seq2)
    R"""
    propspec = 1-($(propGavgc[:,i]))
    propspectrim = propspec[which(propspec > 0)]
    if (length(propspectrim)>0) {
        points(jitter(rep($(seq[seq2[i]]),length(propspectrim))),propspectrim,pch='.',cex=3,col=pal[$i])
    }
    """
end
R"""
lines($(seq[seq2]),1-$mrGavgc,lwd=2)
points($(seq[seq2]),1-$mrGavgc,pch=23,col='black',bg=pal,cex=1.5)
lines($(seq[seq2]),1-$mpGavgc,lwd=2)
points($(seq[seq2]),1 - $mpGavgc, pch=24, col='black', bg=pal, cex=1.5)
#labels
points(4000,1.05,pch=23,col='black',bg='white',cex=1.5,xpd=TRUE)
text(675,1.05,'Functional',cex=1.0,xpd=TRUE)
points(4000,0.95,pch=24,col='black',bg='white',cex=1.5)
text(800,0.95,'Potential',cex=1.0)

#TROPHIC

pal = colorRampPalette(brewer.pal(11,"Spectral"))(length($seq2))
fulldist = $(meantrophicdist[1,:]);
trimdist = fulldist[which(fulldist>0.005)];
plot(trimdist,seq(1,length(trimdist)),type='l',xlim=c(0,0.6),ylim=c(1,12),col=pal[1],xlab='',ylab='',axes=FALSE)
axis(2,at=seq(1:12),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0),las=1)
axis(1,at=seq(0:0.6,by=0.2),labels=TRUE,tck=-0.015,mgp=c(0.5,0.5,0))
title(ylab='Trophic level (TL)', line=2, cex.lab=1.2)
title(xlab='Frequency', line=2.0, cex.lab=1.2)
mtext(paste0("c"), side = 3, adj = -0.3, 
    line = -0.6,cex=1.2,font=2)
points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[1],col='black')
"""
for i=1:12
    R"""
    rect(0,$i-0.5,$(meantrophicdist[length(seq2),i]),$i+0.5,col=paste(pal[length($seq2)],50,sep=''),border=NA)
    """
end
for i=1:length(seq2)
    R"""
    fulldist = $(meantrophicdist[i,:]);
    trimdist = fulldist[which(fulldist>0.005)];
    lines(trimdist,seq(1,length(trimdist)),col=pal[$i])
    points(trimdist,seq(1,length(trimdist)),pch=21,bg=pal[$i],col='black')
    """
end
R"""
legend(x=0.45,y=12.5,legend=$(seq[seq2]),pt.cex=1.2,pt.bg=pal,col='black',pch=21,bty='n',title='',cex=0.8)
text(0.45,11.7,'Assembly time')
dev.off()
"""















