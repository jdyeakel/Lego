function intmatrixv3(S, lambda, probs)
    
    #NOTE In this version, interactions are randomly assigned and there is no inputted distribution for trophic or need interactions
    
    p_n=probs.p_n;
    p_a=probs.p_a;
    # p_i=copy(probs[3]);
    #Defining paiwise probabilities
    #These are probabilities of pairwise interactions within the whole universe of possible interactions (15)

    
    
    #Draw the number of objects made per species
    #Done many times, the mean will be lambda*S
    #A species is an engineer if number of objects > 0
    pdist = Poisson(lambda);
    OpS = rand(pdist,S-1);
    OpS = [0;OpS]; # basal resource does not make anything
    engind = findall(!iszero,OpS);
    E = length(engind);
    
    #Some of these objects may not be unique
    obline = collect(1:sum(OpS));
    obindpS = zeros(Int64,E,sum(OpS));
    for i = 1:E
        o = sample(obline,OpS[engind][i],replace=false);
        obindpS[i,o] .= 1;
    end
    
    obindpS = obindpS[:,findall(!iszero,vec(sum(obindpS,dims=1)))];
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
    
    pwp = (
      na = p_n*(p_a/(p_a+p_n+p_i+p_m)) + p_a*(p_n/(p_a+p_i+p_n)),
      nn = p_n*(p_n/(p_a+p_n+p_i+p_m)),
      ni = p_n*(p_i/(p_a+p_n+p_i+p_m)) + p_i*(p_n/(p_a+p_n+p_i)),
      nm = p_n*(p_m/(p_a+p_n+p_i+p_m)) + p_m*1, #(p_n/p_n),
      ia = p_i*(p_a/(p_a+p_n+p_i)) + p_a*(p_i/(p_a+p_i+p_n)),
      ii = p_i*(p_i/(p_a+p_n+p_i)),
      aa = p_a*(p_a/(p_i+p_n+p_a))
    );

    #Create an empty character array with dimensions equal to the number of players
    #Region 1) Upper left SxS are species-species interactions
    #Region 2) Upper right SxO are species-object interactions
    #Region 3) Lower left OxS are object-species interacions (n only)
    #Region 4) Lower right OxO are object-object interactions (i only)
    int_m = Array{Char}(undef,N,N);
    int_m .= Ref('0');
    
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
    
    #Matrix of trophic-only interactions
    #The number of species in the trophic and mutualistic matrices will all start out the sames size as int_m, however, the objects will be trimmed later
    tp_m = zeros(Int64,S,S);
    #Matrix of mutualistic interactions
    mp_m = zeros(Int64,S,S);

    #The first true species (row/col 2) is always a primary producer
    int_m[2,1] = 'a'
    tp_m[2,1] = 1;
    
    #Fill in diagonal
    #Index 1 is the basal resource ('i')
    #Indices 2:S are species ('i')
    #Indices S+1:N are objects ('i')
    diagindices = diagind(int_m);

    int_m[diagindices[2:S]] .= Ref('n');

    #Deal with the basal resource
    int_m[1,:] .= Ref('i');
    tp_m[1,:] .= 0;
    int_m[findall(x->x=='0',int_m[:,1]),1] .= Ref('i');

    
    #NOTE: Interactions BETWEEN SPECIES

    #Eliminate n-a, i-a, a-a, n-m interactions, which are already determined
    
    # 1) pr_na
    # 2) pr_nn ** 
    # 3) pr_ni **
    # 4) pr_nm
    # 5) pr_ia
    # 6) pr_ii **
    # 7) pr_aa
    # deleteat!(pw_prob_new,[1,4,5,7])
    # pw_prob_new = pw_prob_init[[1,2,3,5,6,7]];
    pw_prob_new = [pwp.na,pwp.nn,pwp.ni,pwp.ia,pwp.ii,pwp.aa];
    pw_prob_new = pw_prob_new/sum(pw_prob_new);


    prob_line = cumsum(pw_prob_new);



    for i = 2:S
        for j = 2:S
            #Only choose interactions for empty elements (int_m = '0')
            #Only do this for the lower triangular part of the matrix
            if int_m[i,j] == '0' && i > j

                r_draw = rand()
                
                #N:A - asymmetric mutualism
                if r_draw < prob_line[1]
                  rr_draw = rand();
                  if rr_draw < 0.5
                      int_m[i,j] = 'a';
                      tp_m[i,j] = 1;
                      int_m[j,i] = 'n';
                      mp_m[j,i] = 1;
                  else
                      int_m[i,j] = 'n';
                      mp_m[i,j] = 1;
                      int_m[j,i] = 'a';
                      tp_m[j,i] = 1;
                  end
                end

                #N:N - symmetric mutualism
                if r_draw > prob_line[1] && r_draw < prob_line[2]
                  int_m[i,j] = 'n'
                  int_m[j,i] = 'n'
                  #Update mutualistic network
                  mp_m[i,j] = 1;
                  mp_m[j,i] = 1;
                end

                #N:I - commensalism
                if r_draw > prob_line[2] && r_draw < prob_line[3]
                    rr_draw = rand();
                    if rr_draw < 0.5
                        int_m[i,j] = 'i';
                        int_m[j,i] = 'n';
                        mp_m[j,i] = 1;
                    else
                        int_m[i,j] = 'n';
                        mp_m[i,j] = 1;
                        int_m[j,i] = 'i';
                    end
                end
                
                #I:A - predation
                if r_draw > prob_line[3] && r_draw < prob_line[4]
                    rr_draw = rand();
                    if rr_draw < 0.5
                        int_m[i,j] = 'a';
                        tp_m[i,j] = 1;
                        int_m[j,i] = 'i';
                    else
                        int_m[i,j] = 'i';
                        int_m[j,i] = 'a';
                        tp_m[j,i] = 1;
                    end
                end
                
                #I:I - neutral interaction
                if r_draw > prob_line[4] && r_draw < prob_line[5]
                  int_m[i,j] = 'i'
                  int_m[j,i] = 'i'
                end
                
                #I:A - predation
                #N:N - symmetric mutualism
                if r_draw > prob_line[5] && r_draw < prob_line[6]
                  int_m[i,j] = 'a'
                  int_m[j,i] = 'a'
                  #Update mutualistic network
                  tp_m[i,j] = 1;
                  tp_m[j,i] = 1;
                end

            end #if
        end #j
    end #i
    
    #We could assume that any species without recorded trophic interactions is a primary producer
    total_trophic = vec(sum(tp_m,dims=2));
    prim_prod = deleteat!(findall(iszero,total_trophic),1); #eliminate row 1
    int_m[prim_prod,1] .= Ref('a');
    tp_m[prim_prod,1] .= 1;
    
    #SPECIES-OBJECT INTERACTIONS
    
    # 1) pr_na
    # 2) pr_nn
    # 3) pr_ni **
    # 4) pr_nm
    # 5) pr_ia
    # 6) pr_ii **
    # 7) pr_aa

    # so_pw_prob = pw_prob_init[[3,5,6]];
    so_pw_prob = [pwp.ni,pwp.ia,pwp.ii];
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
                
                #assimilate-ignore
                if r_draw > so_prob_line[1] && r_draw < so_prob_line[2]
                    int_m[i,j] = 'a';
                    int_m[j,i] = 'i';
                end
                
                #ignore-ignore
                if r_draw > so_prob_line[2] && r_draw < so_prob_line[3]
                    int_m[i,j] = 'i';
                    int_m[j,i] = 'i';
                end
            end
        end
    end
    #All object-object interactions are 'i-i'
    int_m[obind,obind] .= Ref('i');
    
    #A matrix for Direct + Indirect trophic interactions
    tind_m = copy(tp_m);
    #A matrix for Direct + Indirect mutualistic interactions
    mind_m = copy(mp_m);

    #All object-object interactions are 'i-i'
    int_m[obind,obind] .= Ref('i');
    
    #Force Basal primary producer (Row/Col 2) to not 'need' anything
    colonizer_n = deleteat!(findall(x->x=='n',int_m[2,:]),1);
    
    int_m[2,colonizer_n] .= Ref('i');
    mp_m[2,:] .= 0;
    
    #Document the indirect interactions
    for i=2:S
        for j=1:O
            if int_m[spind[i],obind[j]] == 'a'
                makers = findall(x->x=='m',int_m[:,obind[j]]);
                tind_m[spind[i],makers] .= 1;
            end
            if int_m[spind[i],obind[j]] == 'n'
                makers = findall(x->x=='m',int_m[:,obind[j]]);
                mind_m[spind[i],makers] .= 1;
            end
        end
    end
    
    
    return(int_m, tp_m, tind_m, mp_m, mind_m)
    
end
