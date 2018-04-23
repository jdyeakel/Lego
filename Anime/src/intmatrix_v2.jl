function intmatrix(S, probs, ppweight)
    
    p_n=copy(probs[1]);
    p_a=copy(probs[2]);
    p_m=copy(probs[3]);
    p_i=copy(probs[4]); #Ignore with 1 - pr(sum(other))
    #Defining paiwise probabilities
    #These are probabilities of pairwise interactions within the whole universe of possible interactions (15)
    pw_prob_init = [
      pr_na = p_n*(p_a/(p_a+p_n+p_i+p_m)) + p_a*(p_n/(p_a+p_i+p_n)),
      pr_nn = p_n*(p_n/(p_a+p_n+p_i+p_m)),
      pr_ni = p_n*(p_i/(p_a+p_n+p_i+p_m)) + p_i*(p_n/(p_a+p_n+p_i)),
      pr_nm = p_n*(p_m/(p_a+p_n+p_i+p_m)) + p_m*1, #(p_n/p_n),
      pr_ia = p_i*(p_a/(p_a+p_n+p_i)) + p_a*(p_i/(p_a+p_i+p_n)),
      pr_ii = p_i*(p_i/(p_a+p_n+p_i)),
      pr_aa = p_a*(p_a/(p_i+p_n+p_a))
    ];
    
    #Expected size of the system
    O = convert(Int64,round(pr_nm*S*(S-1),0));
    N = S + O;

    #Calculate the distribution of trophic links per species
    #Step 1: niche values
    nichev = rand(S);
    
    #Directed Connectance
    #Connectance for all of the species and objects in the system... so based on num_play and not S
    C = (pr_ia+pr_na+pr_aa); #*(num_play*(num_play-1)))/(num_play^2);
    Ebeta = 2*C;
    beta = (1/Ebeta) - 1;
    BetaD = Beta(1,beta);
    rBeta = rand(BetaD,S);
    rangev = rBeta .* nichev;
    degrees = convert.(Int,round.(rangev * S,0));
    degrees[find(iszero,degrees)]=1;


    #Create an empty character array with dimensions equal to the number of players
    #Region 1) Upper left SxS are species-species interactions
    #Region 2) Upper right SxO are species-object interactions
    #Region 3) Lower left OxS are object-species interacions (n only)
    #Region 4) Lower right OxO are object-object interactions (i only)
    int_m = Array{Char}(N,N);
    #Matrix of trophic-only interactions
    #The number of species in the trophic and mutualistic matrices will all start out the sames size as int_m, however, the objects will be trimmed later
    t_m = zeros(Int64,S,S);
    #Matrix of mutualistic interactions
    m_m = zeros(Int64,S,S);

    prim_prod = zeros(S);
    #The first true species (row/col 2) is always a primary producer
    prim_prod[2] = 1;

    #Assigning trophic interactions to species in the net
    for i = 2:S
        
        #Species trophic degree
        k = degrees[i];
        #Possible foods
        resource = sample(collect(2:N),k,replace=false);
        
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
        t_m[i,spresource] = 1;


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
        m_m[res_n,i] = 1;
        m_m[i,res_n] = 1;

    end

    #Duplicate this information for tall_m
    #This matrix includes indirect trophic interactions
    #Where if 1 (m) 2, and 3 (a) 2, then 3 (a) 1
    tall_m = copy(t_m);


    #Non-trophic interactions

    int_labels = ["na","nn","ni","nm","ia","ii","aa"]

    pw_prob_new = copy(pw_prob_init)
    #Eliminate n-a, i-a, a-a interactions, which are already determined
    deleteat!(pw_prob_new,[1,5,7])
    pw_prob_new = pw_prob_new/sum(pw_prob_new)

    new_int_labels = copy(int_labels)
    deleteat!(new_int_labels,[1,5,7])



    prob_line = cumsum(pw_prob_new);



    for i = 2:num_play
    for j = 2:num_play
      #Only choose interactions for empty elements (int_m = '0')
      #Only do this for the lower triangular part of the matrix
      if int_m[i,j] == '0' && i > j

        r_draw = rand()

        #N:N - symmetric mutualism
        if r_draw < prob_line[1]
          int_m[i,j] = 'n'
          int_m[j,i] = 'n'
          #Update mutualistic network
          m_m[i,j] = 1;
          m_m[j,i] = 1;
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


        #N:M - ecosystem engineering
        if r_draw > prob_line[2] && r_draw < prob_line[3]
          index = 3 #find(x -> x == "nm",new_int_labels)
          mij = Array{Char}(2)
          mij[1] = rand([new_int_labels[index][1],new_int_labels[index][2]])
          if mij[1] == new_int_labels[index][1]
            mij[2] = new_int_labels[index][2]
          else mij[2] = new_int_labels[index][1]
          end
          int_m[i,j] = mij[1]
          int_m[j,i] = mij[2]
          #Ensure nothing 'makes' the forced primary producer
          #Not sure if I need this 10/24/16
          if j == 2
            int_m[i,j] = 'n'
            int_m[j,i] = 'm'
          end
        end


        #I:I - neutral interaction
        if r_draw > prob_line[3] && r_draw < prob_line[4]
          int_m[i,j] = 'i'
          int_m[j,i] = 'i'
        end

      end #if
    end #j
    end #i

    #NOTE: is this necessary??? Could make all 'n' and the objects get turned into 'i's below.
    #(and I think that everything should have at least one 'a')

    #Filling in the diagonal
    #Determine which are active players and which aren't
    for i = 1:num_play
    if length(find(x->x=='a',int_m[i,:])) > 0
      int_m[i,i] = 'n'
    else int_m[i,i] = 'i'
    end
    end


    #Implement posthoc Rules
    #1: Row/Col 1 is the sun
    int_m[1,:] = 'i'
    #All column elements that ARENT 'a' are 'ignore'
    force_ignore=find(x->x!='a',int_m[:,1]);
    int_m[force_ignore,1] = 'i';

    #due to the above forcing, there are no mutualistic interactions with the sun
    #The only 'a' interactions with the sun are 'a-i', not 'a-n'
    m_m[1,:] = 0; #this probably isn't needed.
    m_m[:,1] = 0;

    #A matrix for Direct + Indirect mutualistic interactions
    mall_m = copy(m_m);


    #Basal primary producer (Row/Col 2) doesn't need anything
    colonizer_n=find(x->x=='n',int_m[2,:]);
    int_m[2,colonizer_n[2:length(colonizer_n)]] = 'i';

      