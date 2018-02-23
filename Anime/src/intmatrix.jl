function intmatrix(num_play, probs, ppweight)

    p_n=copy(probs[1]);
    p_a=copy(probs[2]);
    p_m=copy(probs[3]);
    p_i=copy(probs[4]); #Ignore with 1 - pr(sum(other))
    #Defining paiwise probabilities
    #These are probabilities of pairwise interactions within the whole universe of possible interactions (15)
    #THIS IS PROBABLY RIGHT
    pw_prob_init = [
      pr_na = p_n*(p_a/(p_a+p_n+p_i+p_m)) + p_a*(p_n/(p_a+p_i+p_n)),
      pr_nn = p_n*(p_n/(p_a+p_n+p_i+p_m)),
      pr_ni = p_n*(p_i/(p_a+p_n+p_i+p_m)) + p_i*(p_n/(p_a+p_n+p_i)),
      pr_nm = p_n*(p_m/(p_a+p_n+p_i+p_m)) + p_m*1, #(p_n/p_n),
      pr_ia = p_i*(p_a/(p_a+p_n+p_i)) + p_a*(p_i/(p_a+p_i+p_n)),
      pr_ii = p_i*(p_i/(p_a+p_n+p_i)),
      pr_aa = p_a*(p_a/(p_i+p_n+p_a))
    ];

    #Calculate the distribution of trophic links per species


    #Step 1: niche values
    nichev = rand(num_play);
    #Derive connectance from the pr('assimilate')*number of potential links
    #But first try to discount the number of players that are NOT species, which should be on average num_play*pr_nm
    S = num_play - (pr_nm*num_play);
    #Pr('A')*Directed Connectance
    C = ((pr_ia+pr_na+pr_aa)*(S*(S-1)))/(S^2);
    Ebeta = 2*C;
    beta = (1/Ebeta) - 1;
    BetaD = Beta(1,beta);
    rBeta = rand(BetaD,num_play);
    #NOTE: Perhaps I should accurately count how many prey fall within the range
    # rangev = length(find(x->x>nmin&&x<nmax,nichev));
    rangev = rBeta .* nichev;
    centerv = zeros(num_play);
    nmin = zeros(num_play);
    nmax = zeros(num_play);
    degrees = zeros(Int64,num_play);
    for i=1:num_play
    ud = Uniform(rangev[i]/2,nichev[i]);
    centerv[i] = rand(ud);
    nmin[i] = centerv[i] - (rangev[i]/2);
    nmax[i] = centerv[i] + (rangev[i]/2);
    degrees[i] = length(find(x->x>nmin[i]&&x<nmax[i],nichev));
    #Force everything to eat SOMETHING
    if degrees[i] == 0
      degrees[i] = 1;
    end
    end


    #Create an empty character array with dimensions equal to the number of players
    int_m = Array{Char}(num_play,num_play);
    #Matrix of trophic-only interactions
    #The number of species in the trophic and mutualistic matrices will all start out the sames size as int_m, however, the objects will be trimmed later
    t_m = zeros(Int64,num_play,num_play);
    #Matrix of mutualistic interactions
    m_m = zeros(Int64,num_play,num_play);

    #Set array equal to zero
    int_m[1:(num_play*num_play)] = '0';

    prim_prod = zeros(num_play);
    #The first true species (row/col 2) is always a primary producer
    prim_prod[2] = 1;


    #Assigning trophic interactions to species in the net
    for i = 2:num_play

        #Assigning initial trophic interactions
        #We could add an additional vector of 1s to weight towards the basal resource
        #if we use a weighted vector, say 1/3 of species will be primary producers
        pvec = Array{Int64}(Int64(round(num_play*(ppweight))))*0+1;
        vec = collect(1:num_play);
        #Remove ith element
        deleteat!(vec,i);
        #weighted vec
        wvec = [pvec;vec];
        #What is the degree of player i?
        k = degrees[i]
        #Randomly choose the identities of prey without replacement
        if prim_prod[i] == 0
          resource = sample(wvec,k,replace=false)
        else
          #If you are a primary producer, at least 1 thing you conusume must be row/column 1
          resource_wop = sample(wvec,k-1,replace=false)
          resource = [resource_wop;1]
        end


      #Establish these prey in the interaction matrix
      int_m[i,collect(resource)] = 'a'
      #Initially, the resources are all species. Those that get turned into objects after 'm' interactions are deteremined will be trimmed from the trophic matrix (because only interactions with species are allowed)

      #Assigning complimentary interactions
      #First disconsider aa interactions
      aa_int = find(x -> x == 'a', int_m[collect(resource),i])
      aa_int_loc = resource[aa_int];


      if length(aa_int) >= 1
        #Remove elements ee_int
        deleteat!(resource,aa_int)
      end

      #An 'a' has been drawn...
      #Now randomly determine ai or an (aa already emerges from assigning a's)
      #Draw from a binomial: which interaciton will be Asymmetric predation?
      #Draw 1 = Asymmetric predation; Draw 0 = facultative mutualism
      pr_i_given_a = (p_i/(p_i+p_n)); #again discounting a-a, which are already determined
      bindist = Binomial(1,pr_i_given_a)
      bin_draw = rand(bindist,length(resource))

      #Asymmetric predation
      res_i = resource[find(x->x==1,bin_draw)]
      int_m[collect(res_i),i] = 'i' #assigning n's according to random draw

      #Faculatative mutualism
      res_n = resource[find(x->x==0,bin_draw)]
      int_m[collect(res_n),i] = 'n' #assigning i's according to random draw

      #Record mutualistic interaction in m_m matrix
      # m_m[collect(res_n),i] = 1;
      # m_m[i,collect(res_n)] = 1;

    end

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
            # m_m[i,j] = 1;
            # m_m[j,i] = 1;
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


    #Basal primary producer (Row/Col 2) doesn't need anything
    colonizer_n=find(x->x=='n',int_m[2,:]);
    int_m[2,colonizer_n[2:length(colonizer_n)]] = 'i';

    #Assigning objects
    for i = 2:num_play
      if in('m',int_m[:,i])==true
        #The made thing ignores everything
        int_m[i,:] = 'i';
      end
    end
    made = find(x->x=='i',diag(int_m));
    species = find(x->x=='n',diag(int_m));
    for i=2:length(made)
      object = made[i];
      makers = find(x->x=='m',int_m[:,made[i]]);
      #If there is a species that already makes this object, establish m-n
      if length(makers) > 0
        #Reset m-n interactions
        int_m[object,makers] = 'n';
      else
        #Sometimes an object was made by an object, but this gets deleted in the first step. So choose a new species to link the adrift object to.
        #find a species to link it to (i-i to m-n)
        iispecies_all = find(x->x=='i',int_m[:,object]);
        #Take out row 1 and row 'object'... these cannot 'make'
        iispecies = filter(x->x!=1&&x!=object,iispecies_all);
        newmaker = rand(iispecies);
        int_m[object,newmaker] = 'n';
        int_m[newmaker,object] = 'm';
      end


    end

  return(int_m)

  end
