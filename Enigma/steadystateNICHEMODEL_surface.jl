

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


#The range of values to explore
cnrange = collect(0.1:(20-.1)/20:20);
cerange = collect(0.1:(20-.1)/20:20);
# cprange = collect(0.1:(10-.1)/999:10);
# p_nrange = collect(0.0005:(0.003-0.0005)/999:0.003);
# p_arange = collect(0.005:(0.03-0.005)/999:0.03);

error = Array{Float64}(undef,size(cnrange)[1],size(cerange)[1]);

# cn = cnvec[208];
# ce = cevec[134];
# cp = cpvec[92];
# p_n = p_nvec[11];
# p_a = p_avec[75];
# errmean_old = 5.;
# temperature = 1.0;
# tic = 1;

cn_init = 3.14;
ce_init = 1.41;
cp_init = 1.0;
p_n_init = 0.002;
p_a_init = 0.01;

cp = 1;
p_n = p_n_init;
p_a = p_a_init;
errmean_old = 1.;
temperature = 1.0;
tic = 1;

    # global cn = 4.195216435170081;
    # global ce = 1.8682437375935583;
    # global cp = 0.6511860914279262;
    # global p_n = 0.0030274490082625636;
    # global p_a =  0.014690692329311795;


for r=1:size(cnrange)[1]
    for s=1:size(cerange)[1]

        # tempvec[r] = copy(temperature);
        
        #Select proposed value from predefined range with uniform selection over a smaller range determined by temperature
        # cn_prop = proposal(temperature/1,cn,cnrange);
        # ce_prop = proposal(temperature/1,ce,cerange);
        # cp_prop = proposal(temperature/1,cp,cprange);
        # p_n_prop = proposal(temperature/2,p_n,p_nrange);
        # p_a_prop = proposal(temperature/2,p_a,p_arange);
        
        cn = cnrange[r];
        ce = cerange[s];
        
        S = 200;
        maxits = 2000;
        SOprobs = (
        p_n=p_n_init,
        p_a=p_a_init
        );
        SSmult = 1.0; OOmult = 0.0;
        SSprobs = (p_n = SSmult .* SOprobs.p_n , p_a = SSmult .* SOprobs.p_a);
        OOprobs = (p_n = OOmult * SOprobs.p_n, p0 = 0.0);


        #expected objects per species
        lambda = 0;
        athresh = 0;
        nthresh = 1.0;
        MaxN = convert(Int64,floor(S + S*lambda));

        enigmareps=50;
        jspecies = SharedArray{Float64}(enigmareps);
        jconn = SharedArray{Float64}(enigmareps);
        jmdegree = SharedArray{Float64}(enigmareps);
        jstdindegree = SharedArray{Float64}(enigmareps);
        jstdoutdegree = SharedArray{Float64}(enigmareps);
        @sync @distributed for j = 1:enigmareps
            int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv4(S,lambda,SSprobs,SOprobs,OOprobs);

            a_b,
            n_b,
            i_b,
            m_b,
            n_b0,
            sp_v,
            int_id = preamble_defs(int_m);

            sprich,rich,clock,CID = assembly(
                int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
                athresh,nthresh,maxits,cn_prop,ce_prop,cp_prop);

            Aenigma = a_b[CID[:,maxits],CID[:,maxits]];
            Senigma = size(Aenigma)[1];
            Aenigma = Aenigma[2:Senigma,2:Senigma];
            Senigma -= 1;
            Cenigma = sum(Aenigma)/(Senigma^2);

            jspecies[j] = Senigma;
            jconn[j] = Cenigma;
            jmdegree[j] = mean(sum(Aenigma,dims=2));
            jstdindegree[j] = std(sum(Aenigma,dims=2));
            jstdoutdegree[j] = std(sum(Aenigma,dims=1));
        end



        enig_species = mean(jspecies[vec(findall(!isnan,jspecies))]);
        enig_conn = mean(jconn[vec(findall(!isnan,jconn))]);
        enig_mdegree = mean(jmdegree[vec(findall(!isnan,jmdegree))]);
        enig_stdindegree = mean(jstdindegree[vec(findall(!isnan,jstdindegree))]);
        enig_stdoutdegree = mean(jstdoutdegree[vec(findall(!isnan,jstdoutdegree))]);


        #Simulate a bunch of nichemodelwebs



        nichereps = 100;
        ispecies = Array{Float64}(undef,nichereps);
        iconn = Array{Float64}(undef,nichereps);
        imdegree = Array{Float64}(undef,nichereps);
        istdindegree = Array{Float64}(undef,nichereps);
        istdoutdegree = Array{Float64}(undef,nichereps);

        for i=1:nichereps
            #NICHE MODEL COMPARISON
            #NOTE: C CANNOT BE EQUAL TO OR GREATER THAN 0.5
            #NOTE Violates beta parameters!
            Aniche, n = nichemodelweb(Int64(floor(species[r,2])),conn[r,2]);

            #make measurements
            Aniche = Array{Int64}(Aniche);

            Sniche = size(Aniche)[1];
            Cniche = sum(Aniche)/Sniche^2;

            ispecies[i] = Sniche;
            #Connectance
            iconn[i] = Cniche;

            #mean degree
            imdegree[i] = mean(sum(Aniche,dims=2));

            #indegree
            istdindegree[i] = std(sum(Aniche,dims=2));

            #outdegree
            istdoutdegree[i] = std(sum(Aniche,dims=1));
        end
        # 
        #ignore NAN
        
        niche_species = mean(ispecies);
        
        #Connectance
        niche_conn = mean(iconn);
        
        #mean degree
        niche_mdegree = mean(imdegree);
        
        
        #indegree
        niche_stdindegree = mean(istdindegree);
        
        
        #outdegree
        niche_stdoutdegree = mean(istdoutdegree);


        #Calculate error
        # err_sp = sqrt((species[r,2] - species[r,1])^2)/mean(species[r,1]); #/std(ispecies);
        # err_conn = sqrt((conn[r,2] - conn[r,1])^2)/mean(conn[r,1]); #/std(iconn);
        # err_md = sqrt((mdegree[r,2] -mdegree[r,1])^2)/mean(mdegree[r,1]); #/std(imdegree);
        # err_sdin = sqrt((stdindegree[r,2] - stdindegree[r,1])^2)/mean(stdindegree[r,1]); #/std(istdindegree);
        # err_sdout = sqrt((stdoutdegree[r,2] - stdoutdegree[r,1])^2)/mean(stdoutdegree[r,1]); #/std(istdoutdegree);

        error[r,s] = sqrt((enig_mdegree - niche_mdegree)^2); #/std(imdegree);
        println(string("r=",r,"; s=",s))
    end
end


filename = "figures/niche/surface.pdf"
namespace = smartpath(filename);
R"""
library(RColorBrewer)
pal = brewer.pal(5,'Set1')
pdf($namespace,height=5,width=6)
image(x=$cnrange,y=$cerange,z=$(error))
"""
R"dev.off()"
