# loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");
loadfunc = include("$(homedir())/2014_Lego/Anime/src/loadfuncsYOG.jl");


S = 400;

tmax = 2000;
tswitch = 1000;

# S = 400;
probs = [
p_n=0.003,
p_a=0.003
# p_n = 0.02,
# p_a = 0.02
];
#expected objects per species
lambda = 0.5;

a_thresh = 0;
n_thresh = 0.2;

extinctions = [ones(Float64,tswitch);ones(Float64,tmax-tswitch)];
colonizations = [ones(Float64,tswitch);ones(Float64,tmax-tswitch)];

lvec = 10;
#Predator load parameters
epsilonvec = collect(0.01:(0.5-0.01)/(lvec-1):0.5);
sigmavec = collect(0.1:(0.5-0.1)/(lvec-1):0.5);

#Resource overlap parameters
extmidvec = collect(0.1:(1.0-0.1)/(lvec-1):1.0);
steepvec = collect(0.5:(3.0-0.5)/(lvec-1):3.0);

reps = 50;

ROSS = Array{Float64}(length(extmidvec),length(steepvec));

PLSS = Array{Float64}(length(epsilonvec),length(sigmavec));
#Search across the parameter space for both PL and RO models and take average steady state across repetitions
for i = 1:length(extmidvec)
    for j = 1:length(steepvec)
        
        extmid = extmidvec[i];
        steep = steepvec[j];
        
        epsilon = epsilonvec[i];
        sigma = sigmavec[j];
        
        rss = SharedArray{Float64}(reps);
        pss = SharedArray{Float64}(reps);
        
        @sync @parallel for r = 1:reps

            int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

            a_b,
            n_b,
            i_b,
            m_b,
            n_b0,
            sp_v,
            int_id = preamble_defs(int_m);
            
            exttype = "RO"; #RO #PL #ROPL
            
            sprich,
            cid,
            lpot_col,
            status,
            prim_ext,
            sec_ext,
            CID = assembly_trim(
                int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
                a_thresh,n_thresh,extmid,steep,epsilon,sigma,colonizations,extinctions,tmax,exttype);
            
            avgss = mean(sprich[tmax-100:tmax]);
            rss[r] = avgss;
                
            exttype = "PL"; #RO #PL #ROPL
            
            sprich,
            cid,
            lpot_col,
            status,
            prim_ext,
            sec_ext,
            CID = assembly_trim(
                int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
                a_thresh,n_thresh,extmid,steep,epsilon,sigma,colonizations,extinctions,tmax,exttype);
            
            avgss = mean(sprich[tmax-100:tmax]);
            pss[r] = avgss;
            
        end
        
        ROSS[i,j] = mean(rss);
        PLSS[i,j] = mean(pss);
        
    end
end

            

save(string("$(homedir())/2014_Lego/Anime/data/SS.jld"),
"ROSS",ROSS,
"PLSS",PLSS,
"epsilonvec",epsilonvec,
"sigmavec",sigmavec,
"extmidvec",extmidvec,
"steepvec",steepvec
);

# 
# d = load(string("$(homedir())/2014_Lego/Anime/data/SS.jld"));
# ROSS = d["ROSS"];
# PLSS = d["PLSS"];
# epsilonvec = d["epsilonvec"];
# sigmavec = d["sigmavec"];
# extmidvec = d["extmidvec"];
# steepvec = d["steepvec"];
# 
# 
# namespace = string("$(homedir())/2014_Lego/Anime/figures2/nodegreedist/SS.pdf");
# R"""
# pdf($namespace,height=8,width=18)
# par(mfrow=c(1,2))
# image(x=$extmidvec,y=$steepvec,z=$ROSS)
# image(x=$epsilonvec,y=$sigmavec,z=$PLSS)
# dev.off()
# """
