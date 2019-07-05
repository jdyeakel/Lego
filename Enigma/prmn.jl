if homedir() == "/home/z840"
    loadfunc = include("$(homedir())/2014_Lego/Enigma/src/loadfuncs.jl");
else
    loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Enigma/src/loadfuncs.jl");
end
#git

lambdavec = collect(0:0.05:4.0);
llamb = length(lambdavec);


S = 200;

maxits = 4000;

cn = pi;
ce = sqrt(2);
cp = 1.0;

SOprobs = (
p_n=0.002,
p_a=0.01
);
SSmult = 1.0; OOmult = 0.0;
SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


#expected objects per species
# lambda = 0.0;
athresh = 0;
nthresh = 1.0;
# MaxN = convert(Int64,floor(S + S*lambda));

reps = 1000;

numobj = SharedArray{Float64}(llamb,reps);
numsharedobj = SharedArray{Float64}(llamb,reps);
probsp = SharedArray{Float64}(llamb,reps);

its = llamb*reps;
@sync @distributed for i = 0:(its - 1)

    #Across lambdavec
    a = Int64(floor(i/reps)) + 1;
    #Across reps
    b = mod(i,reps) + 1;

    lambda = lambdavec[a];

    
    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);
            
    #what is the proportion of objects that have multiple makers?
    numobj[a,b] = size(int_m)[1] - S - 1;
    numsharedobj[a,b] = sum(sum(m_b[:,S+1:size(int_m)[1]],dims=1) .> 1);
    
    #If we select 2 species at random, what is the probability that they make the same object?
    cits = 10000;
    let ct = 0
        for i=1:cits
            twospecies = rand(sp_v,2);
            numshare = sum(sum(m_b[twospecies,S+1:size(int_m)[1]],dims=1) .> 1)
            if numshare > 0
                ct += 1;
            end
        end
        probsp[a,b] = ct/cits;
    end
    
end

numuniqueobj = numobj .- numsharedobj;
probsharedobj =  numsharedobj ./ numobj;

m_propsharedobj = mean(probsharedobj,dims=2);
m_numobjects = mean(numobj,dims=2);
m_numsharedobj = mean(numsharedobj,dims=2);
m_numuniqueobj = mean(numuniqueobj,dims=2);

m_probsp = mean(probsp,dims=2);

filename = "data/prmn.jld";
namespace = smartpath(filename);
@save namespace lambdavec llamb reps numobj numsharedobj numuniqueobj probsharedobj probsp;

filename = "data/prmn.jld";
namespace = smartpath(filename);
@load namespace lambdavec llamb reps numobj numsharedobj numuniqueobj probsharedobj probsp;



filename = "../manuscript/fig_redundancy.pdf";
namespace = smartpath(filename);
R"""
pdf($namespace,width=8,height=8)
par(mfrow=c(2,1),oma=c(3,3,1,1),mar=c(1,1,0,1))
boxplot($(transpose(numsharedobj)),names=$lambdavec,xlab='',ylab='Redundant objects',pch='.',col='gray',xaxt='n',yaxt='n')
axis(2)
title(ylab='Number redundant objects',outer=FALSE,xpd=NA,line=2.5)
boxplot($(transpose(numsharedobj ./ numobj)),names=$lambdavec,xlab=expression(paste('Number objects/species (',eta,')')),ylab='',pch='.',col='gray',ylim=c(0,1),xaxt='n',yaxt='n')
axis(1,at=seq(1,length($lambdavec),length.out=5),labels=seq(0,2,by=0.5))
axis(2)
title(xlab=expression(paste('Number objects/species (',eta,')')),outer=TRUE,xpd=NA,line=1.5)
title(ylab='Proportion redundant objects',outer=FALSE,xpd=NA,line=2.5)
dev.off()
"""

