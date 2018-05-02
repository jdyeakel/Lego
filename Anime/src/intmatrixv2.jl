function intmatrixv2(S, lambda, probs, ppweight)
    
    p_n=copy(probs[1]);
    p_a=copy(probs[2]);
    # p_i=copy(probs[3]);
    #Defining paiwise probabilities
    #These are probabilities of pairwise interactions within the whole universe of possible interactions (15)

    
    
    #Draw the number of objects made per species
    #Done many times, the mean will be lambda*S
    #A species is an engineer if number of objects > 0
    pdist = Poisson(lambda);
    OpS = rand(pdist,S-1);
    OpS = [0;OpS]; # basal resource does not make anything
    engind = find(!iszero,OpS);
    E = length(engind);
    
    #Some of these objects may not be unique
    obline = collect(1:sum(OpS));
    obindpS = zeros(Int64,E,sum(OpS));
    for i = 1:E
        o = sample(obline,OpS[engind][i],replace=false);
        obindpS[i,o] = 1;
    end
    
    obindpS = obindpS[:,find(!iszero,sum(obindpS,1))];
    O = size(obindpS)[2];
    
    #Expected size of the system
    # O = convert(Int64,round(lambda*S,0));
    N = S + O;
    spind = collect(1:S);
    obind = collect(S+1:N);
    engind = collect(S-E+1:S);
    
    p_m = (sum(obindpS))/(N^2);
    p_i = 1 - (p_n + p_a + p_m);
    
    
    # exp_p_m = S*lambda/((S + S*lambda*(1-exp(-1)))^2)
    
    pw_prob_init = [
      pr_na = p_n*(p_a/(p_a+p_n+p_i+p_m)) + p_a*(p_n/(p_a+p_i+p_n)),
      pr_nn = p_n*(p_n/(p_a+p_n+p_i+p_m)),
      pr_ni = p_n*(p_i/(p_a+p_n+p_i+p_m)) + p_i*(p_n/(p_a+p_n+p_i)),
      pr_nm = p_n*(p_m/(p_a+p_n+p_i+p_m)) + p_m*1, #(p_n/p_n),
      pr_ia = p_i*(p_a/(p_a+p_n+p_i)) + p_a*(p_i/(p_a+p_i+p_n)),
      pr_ii = p_i*(p_i/(p_a+p_n+p_i)),
      pr_aa = p_a*(p_a/(p_i+p_n+p_a))
    ];

    #Create an empty character array with dimensions equal to the number of players
    #Region 1) Upper left SxS are species-species interactions
    #Region 2) Upper right SxO are species-object interactions
    #Region 3) Lower left OxS are object-species interacions (n only)
    #Region 4) Lower right OxO are object-object interactions (i only)
    int_m = Array{Char}(N,N);
    int_m[:] = '0';
    
    #Assign the make-need interactions
    #Maybe faster way to do this?
    for i = 1:E
        for j=1:O
            if obindpS[i,j] == 1
                int_m[engind[i],obind[j]] = 'm';
                int_m[obind[j],engind[i]] = 'n';
            end
        end
    end
    
    #Calculate the distribution of trophic links per species
    
    #Directed Connectance
    #Connectance for all of the species and objects in the system... so based on num_play and not S
    C = (pr_ia+pr_na+pr_aa); #*(num_play*(num_play-1)))/(num_play^2);
    Ebeta = 2*C;
    beta = (1/Ebeta) - 1;
    BetaD = Beta(1,beta);
    rBeta = rand(BetaD,S);
    nichev = rand(S);
    rangev = rBeta .* nichev;
    degrees = convert.(Int,round.(rangev * S,0));
    degrees[find(iszero,degrees)]=1;
    
    #Matrix of trophic-only interactions
    #The number of species in the trophic and mutualistic matrices will all start out the sames size as int_m, however, the objects will be trimmed later
    tp_m = zeros(Int64,S,S);
    #Matrix of mutualistic interactions
    mp_m = zeros(Int64,S,S);

    prim_prod = zeros(S);
    #The first true species (row/col 2) is always a primary producer
    prim_prod[2] = 1;
    
    #Fill in diagonal
    #Index 1 is the basal resource ('i')
    #Indices 2:S are species ('i')
    #Indices S+1:N are objects ('i')
    diagindices = diagind(int_m);

    int_m[diagindices[2:S]] = 'n';

    #Assigning trophic interactions to species in the net
    for i = 2:S
        
        #Species trophic degree
        k = degrees[i];
        #Possible foods (do not count self or objects you make)
        choices = collect(1:N);
        eliminate = [i;find(x->x=='m',int_m[i,:])];
        deleteat!(choices,eliminate);
        
        resource = sample(choices,k,replace=false);
        
        #Is this species a primary producer? If so, make it eat basal resource
        #Ensure species 2 is always a primary producer
        pp_draw = rand();
        if i != 2
            if pp_draw < ppweight
                resource[rand(collect(1:k))] = 1;
            end
        else 
            resource[rand(collect(1:k))] = 1;
        end
        
        #Which of these resources are species?
        spresource = resource[resource.<=S];
        #Which of these resources are objects?
        obresource = setdiff(resource,spresource);

        #Establish these prey in the interaction matrix
        int_m[i,resource] = 'a'
        tp_m[i,spresource] = 1;


        #Note: to vectorize a row from int_m (so that it is an Array{Char,1}), we would write v = int_m[i,:][:]

        #Assigning complimentary interactions
        #ONLY for interactions with SPECIES
        #First disconsider aa interactions
        aa_int = find(x -> x == 'a', int_m[spresource,i]);
        aa_int_loc = spresource[aa_int];


        if length(aa_int) >= 1
          #Remove elements aa_int
          deleteat!(spresource,aa_int);
        end

        #An 'a' has been drawn...
        #Now randomly determine ai or an (aa already emerges from assigning a's)
        #Draw from a binomial: which interaciton will be Asymmetric predation?
        #Draw 1 = Asymmetric predation; Draw 0 = facultative mutualism
        pr_i_given_a = (p_i/(p_i+p_n)); #again discounting a-a, which are already determined
        bindist = Binomial(1,pr_i_given_a);
        bin_draw = rand(bindist,length(spresource));
        
        #Force ignore for primaru producer basal resource interaction
        ignorebasal = find(x->x==1,spresource);
        if length(ignorebasal) > 0
            bin_draw[ignorebasal] = 1;
        end

        #Asymmetric predation
        res_i = spresource[find(x->x==1,bin_draw)];
        int_m[res_i,i] = 'i'; #assigning i's according to random draw

        #All complimentary interactions for species-object assimilations are 'i'
        int_m[obresource,i] = 'i';

        #Faculatative mutualism
        res_n = spresource[find(x->x==0,bin_draw)];
        int_m[res_n,i] = 'n'; #assigning n's according to random draw

        #Record mutualistic interaction in m_m matrix
        mp_m[res_n,i] = 1;
        # m_m[i,res_n] = 1;

    end

    #Duplicate this information for tind_m
    #This matrix includes indirect trophic interactions
    #Where if 1 (m) 2, and 3 (a) 2, then 3 (a) 1
    tind_m = copy(tp_m);
    
    #Deal with the basal resource
    int_m[1,:] = 'i';
    int_m[find(x->x=='0',int_m[:,1]),1] = 'i';

    
    #NOTE: Non-trophic interactions BETWEEN SPECIES

    int_labels = ["na","nn","ni","nm","ia","ii","aa"]

    pw_prob_new = copy(pw_prob_init)
    #Eliminate n-a, i-a, a-a, n-m interactions, which are already determined
    
    # 1) pr_na
    # 2) pr_nn ** 
    # 3) pr_ni **
    # 4) pr_nm
    # 5) pr_ia
    # 6) pr_ii **
    # 7) pr_aa
    deleteat!(pw_prob_new,[1,4,5,7])
    pw_prob_new = pw_prob_new/sum(pw_prob_new)

    new_int_labels = copy(int_labels)
    deleteat!(new_int_labels,[1,4,5,7])



    prob_line = cumsum(pw_prob_new);



    for i = 2:S
        for j = 2:S
            #Only choose interactions for empty elements (int_m = '0')
            #Only do this for the lower triangular part of the matrix
            if int_m[i,j] == '0' && i > j

                r_draw = rand()

                #N:N - symmetric mutualism
                if r_draw < prob_line[1]
                  int_m[i,j] = 'n'
                  int_m[j,i] = 'n'
                  #Update mutualistic network
                  mp_m[i,j] = 1;
                  mp_m[j,i] = 1;
                end

                #N:I - commensalism
                if r_draw > prob_line[1] && r_draw < prob_line[2]
                  index = 2 #find(x -> x == "ni",new_int_labels)
                  mij = Array{Char}(2)
                  mij[1] = rand([new_int_labels[index][1],new_int_labels[index][2]])
                  if mij[1] == new_int_labels[index][1]
                    mij[2] = new_int_labels[index][2]
                  else mij[2] = new_int_labels[index][1]
                  end
                  int_m[i,j] = mij[1]
                  int_m[j,i] = mij[2]
                end

                #I:I - neutral interaction
                if r_draw > prob_line[2] && r_draw < prob_line[3]
                  int_m[i,j] = 'i'
                  int_m[j,i] = 'i'
                end

            end #if
        end #j
    end #i
    
    # tp_m = copy(t_m);
    # mp_m = copy(m_m);

    #A matrix for Direct + Indirect mutualistic interactions
    mind_m = copy(mp_m);

    
    #Species-object & object-species interactions
    #Draw 'make interaction' probability for all species-objects
    #Assign opposing 'need' interaction
    #Make sure all species produce objects
    
    #Types of non-trophic interactions between species-objects: (a-i already filled in)
    #1) make-need
    #2) need-ignore
    #3) ignore-ignore
    
    #every object needs an enginner
    #This will ensure that each object has a maker... more makers can then be created
    
    
    #This biases the n-m draw by the number of Objects
    #We need to correct for this
    #how many n-m interactions are expected in SxO = pr_nm*S*O
    
    # 1) pr_na
    # 2) pr_nn
    # 3) pr_ni **
    # 4) pr_nm
    # 5) pr_ia
    # 6) pr_ii **
    # 7) pr_aa
    
    so_pw_prob = [pw_prob_init[3],pw_prob_init[6]];
    so_pw_prob = so_pw_prob/sum(so_pw_prob);
    so_prob_line = cumsum(so_pw_prob);
    
    for i=2:S
        for j=S+1:N
            if int_m[i,j] == '0'
                r_draw = rand();
                
                #need-ignore
                if r_draw < so_prob_line[1]
                    int_m[i,j] = 'n';
                    int_m[j,i] = 'i';
                end
                
                # #make-need
                # if r_draw > so_prob_line[1] && r_draw < so_prob_line[2] 
                #     int_m[i,j] = 'm';
                #     int_m[j,i] = 'n';
                # end
                
                #ignore-ignore
                if r_draw > prob_line[1]
                    int_m[i,j] = 'i';
                    int_m[j,i] = 'i';
                end
            end
        end
    end
    
    #All object-object interactions are 'i-i'
    int_m[obind,obind] = 'i';
    
    #Force Basal primary producer (Row/Col 2) to not 'need' anything
    colonizer_n=deleteat!(find(x->x=='n',int_m[2,:]),1);
    int_m[2,colonizer_n] = 'i';
    
    #Document the indirect interactions
    for i=2:S
        for j=1:O
            if int_m[spind[i],obind[j]] == 'a'
                makers = find(x->x=='m',int_m[:,obind[j]]);
                tind_m[spind[i],makers] = 1;
            end
            if int_m[spind[i],obind[j]] == 'n'
                makers = find(x->x=='m',int_m[:,obind[j]]);
                mind_m[spind[i],makers] = 1;
            end
        end
    end
    
    
    return(int_m, tp_m, tind_m, mp_m, mind_m)
    
end
