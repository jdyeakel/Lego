function intmatrixv4(S, lambda, SSprobs, SOprobs, OOprobs)
    
		#NOTE: In this version, different sets of probabilities are used to fill up interactions in each matrix quadrant (SS, SO, OO)
		
    #NOTE In this version, interactions are randomly assigned and there is no inputted distribution for trophic or need interactions
    
		#The SO quadrant has all e-n-i-m interactions, so do this one first
		
    SOp_n=SOprobs.p_n;
    SOp_a=SOprobs.p_a;
    
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
    
    SOp_m = (sum(obindpS))/(N^2);
    SOp_i = 1 - (SOp_n + SOp_a + SOp_m);
    
    
    # exp_p_m = S*lambda/((S + S*lambda*(1-exp(-1)))^2)
    
		#SO interactions: E-N-I-M
		
    SOpwp = (
      SOna = SOp_n*(SOp_a/(SOp_a+SOp_n+SOp_i+SOp_m)) + SOp_a*(SOp_n/(SOp_a+SOp_i+SOp_n)),
      SOnn = SOp_n*(SOp_n/(SOp_a+SOp_n+SOp_i+SOp_m)),
      SOni = SOp_n*(SOp_i/(SOp_a+SOp_n+SOp_i+SOp_m)) + SOp_i*(SOp_n/(SOp_a+SOp_n+SOp_i)),
      SOnm = SOp_n*(SOp_m/(SOp_a+SOp_n+SOp_i+SOp_m)) + SOp_m*1, #(SOp_n/SOp_n),
      SOia = SOp_i*(SOp_a/(SOp_a+SOp_n+SOp_i)) + SOp_a*(SOp_i/(SOp_a+SOp_i+SOp_n)),
      SOii = SOp_i*(SOp_i/(SOp_a+SOp_n+SOp_i)),
      SOaa = SOp_a*(SOp_a/(SOp_i+SOp_n+SOp_a))
    );
		
		#SS interactions: E-N-I
		SSp_n = SSprobs.p_n;
    SSp_a = SSprobs.p_a;
		SSp_i = 1 - SSp_a - SSp_n;
		
		SSpwp = (
      SSna = SSp_n*(SSp_a/(SSp_a+SSp_n+SSp_i)) + SSp_a*(SSp_n/(SSp_a+SSp_i+SSp_n)),
      SSnn = SSp_n*(SSp_n/(SSp_a+SSp_n+SSp_i)),
      SSni = SSp_n*(SSp_i/(SSp_a+SSp_n+SSp_i)) + SSp_i*(SSp_n/(SSp_a+SSp_n+SSp_i)),
      SSia = SSp_i*(SSp_a/(SSp_a+SSp_n+SSp_i)) + SSp_a*(SSp_i/(SSp_a+SSp_i+SSp_n)),
      SSii = SSp_i*(SSp_i/(SSp_a+SSp_n+SSp_i)),
      SSaa = SSp_a*(SSp_a/(SSp_i+SSp_n+SSp_a))
    );
		
		#OO interactions: N-I
		OOp_n = OOprobs.p_n;
		OOp_i = 1 - OOp_n;
		
		OOpwp = (
      OOnn = OOp_n*(OOp_n/(OOp_n+OOp_i)),
      OOni = OOp_n*(OOp_i/(OOp_n+OOp_i)) + OOp_i*(OOp_n/(OOp_n+OOp_i)),
      OOii = OOp_i*(OOp_i/(OOp_n+OOp_i))
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
    pw_prob_new = [SSpwp.na,SSpwp.nn,SSpwp.ni,SSpwp.ia,SSpwp.ii,SSpwp.aa];
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
    so_pw_prob = [SOpwp.ni,SOpwp.ia,SOpwp.ii];
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
		
		    #Object-OBJECT INTERACTIONS
    
    oo_pw_prob = [OOpwp.ni,OOpwp.nn,OOpwp.ii];
    oo_pw_prob = oo_pw_prob/sum(oo_pw_prob);
    oo_prob_line = cumsum(oo_pw_prob);
    
    for i=S+1:N
        for j=S+1:N
            if int_m[i,j] == '0'
                r_draw = rand();
                
                #need-ignore
                if r_draw < oo_prob_line[1]
                    int_m[i,j] = 'n';
                    int_m[j,i] = 'i';
                end
                
                #need-need
                if r_draw > oo_prob_line[1] && r_draw < oo_prob_line[2]
                    int_m[i,j] = 'n';
                    int_m[j,i] = 'n';
                end
                
                #ignore-ignore
                if r_draw > oo_prob_line[2] && r_draw < oo_prob_line[3]
                    int_m[i,j] = 'i';
                    int_m[j,i] = 'i';
                end
            end
        end
    end
    
    #A matrix for Direct + Indirect trophic interactions
    tind_m = copy(tp_m);
    #A matrix for Direct + Indirect mutualistic interactions
    mind_m = copy(mp_m);
    
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
