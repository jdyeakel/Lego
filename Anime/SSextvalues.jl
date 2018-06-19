loadfunc = include("$(homedir())/Dropbox/PostDoc/2014_Lego/Anime/src/loadfuncs.jl");

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

#Predator load parameters
epsilonvec = collect(0.01:0.01:0.1);
sigmavec = collect(0.1:0.1:1.0);

#Resource overlap parameters
extmidvec = collect(0.1:0.1:1.0);
steepvec = collect(0.5:0.16:2.0);

reps = 100;

ROSS = Array{Float64}(length(extmidvec),length(steepvec));

PLSS = Array{Float64}(length(epsilonvec),length(sigmavec));
#Search across the parameter space for both PL and RO models and take average steady state across repetitions
for i = 1:length(extmidvec)
    for j = 1:length(steepvec)
        
        extmid = extmidvec[i];
        steep = steepvec[j];
        
        epsilon = epsilonvec[i];
        sigma = sigmavec[j];
        
        rss = Array{Float64}(reps);
        pss = Array{Float64}(reps);
        
        for r = 1:reps

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

            

            
        
        
