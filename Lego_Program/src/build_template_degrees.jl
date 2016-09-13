function build_template_degrees(num_play, probs)


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
  #pw_prob_init[find(x->x==true,isnan(pw_prob_init))]=0;


  # pw_prob = copy(pw_prob_init) / sum(pw_prob_init)


  # Distributing trophic interactions (e's) according to a degree distribution
  #probabilities are for the whole matrix
  #this means we should have N.e<-p.e*(num_play*(num_play-1)) trophic interactions (e)
  #thus the mean degree should be mean.k=N.e/num_play
  #we can sample degrre from a exponential ditribution with mean (1/rate) equal to mean.k; degrees=rexp(num_play,1/mean.k)


  #Proposed method:
  #Build a random food web from the niche model given the number of species
  #Determine link-density properties from the web
  #Populate the interaction matrix with these statistical properties, but the structure will be scrambled. (I think)


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
  degrees = Array{Int64}(num_play)*0;
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

  #The number of prey in the range should just scale with the number of 'species'
  #degreesAn = rangev .* num_play;
  #The +1 forces species to eat at least one thing
  #degreesAn = round(Int64,copy(degreesAn)) + 1;
  #plot(x=degrees,Geom.histogram)

  #Create an empty character array with dimensions equal to the number of players
  int_m = Array(Char,num_play,num_play);
  t_mpre = zeros(num_play,num_play);
  #Set array equal to zero
  int_m[1:(num_play*num_play)] = '0';

  prim_prod = zeros(num_play);
  prim_prod[sample(collect(1:num_play),1)] = 1;

  for i = 2:num_play

    #Assigning initial trophic interactions
    #We could add an additional vector of 1s to weight towards the basal resource
    vec = collect(1:num_play);
    #Remove ith element
    deleteat!(vec,i);
    #What is the degree of player i?
    k = degrees[i]
    #Randomly choose the identities of prey without replacement
    if prim_prod[i] == 0
      resource = sample(vec,k,replace=false)
    else
      resource_wop = sample(vec,k-1,replace=false)
      resource = [resource_wop;1]
    end


    #Establish these prey in the interaction matrix
    int_m[i,collect(resource)] = 'a'
    t_mpre[i,collect(resource)] = 1;
    t_mpre[collect(resource),i] = 1;
    #Note: to vectorize a row from int_m (so that it is an Array{Char,1}), we would write v = int_m[i,:][:]

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

  end


  #Non-trophic interactions

  int_labels = ["na","nn","ni","nm","ia","ii","aa"]

  pw_prob_new = copy(pw_prob_init)
  #Eliminate n-a, i-a, a-a interactions, which are already determined
  deleteat!(pw_prob_new,[1,5,7])
  pw_prob_new = pw_prob_new/sum(pw_prob_new)

  new_int_labels = copy(int_labels)
  deleteat!(new_int_labels,[1,5,7])



  prob_line = cumsum(pw_prob_new) #Why did I sort this???



  for i = 2:num_play
    for j = 2:num_play
      #Only choose interactions for empty elements (int_m = '0')
      #Only do this for the lower triangular part of the matrix
      if int_m[i,j] == '0' && i > j

        r_draw = rand()

        #N:N
        if r_draw < prob_line[1]
          int_m[i,j] = 'n'
          int_m[j,i] = 'n'
        end

        #N:I
        if r_draw > prob_line[1] && r_draw < prob_line[2]
          index = 2 #find(x -> x == "ni",new_int_labels)
          mij = Array(Char,2)
          mij[1] = rand([new_int_labels[index][1],new_int_labels[index][2]])
          if mij[1] == new_int_labels[index][1]
            mij[2] = new_int_labels[index][2]
          else mij[2] = new_int_labels[index][1]
          end
          int_m[i,j] = mij[1]
          int_m[j,i] = mij[2]
        end

        #N:M
        if r_draw > prob_line[2] && r_draw < prob_line[3]
          index = 3 #find(x -> x == "nm",new_int_labels)
          mij = Array(Char,2)
          mij[1] = rand([new_int_labels[index][1],new_int_labels[index][2]])
          if mij[1] == new_int_labels[index][1]
            mij[2] = new_int_labels[index][2]
          else mij[2] = new_int_labels[index][1]
          end
          int_m[i,j] = mij[1]
          int_m[j,i] = mij[2]
        end


        #I:I
        if r_draw > prob_line[3] && r_draw < prob_line[4]
          int_m[i,j] = 'i'
          int_m[j,i] = 'i'
        end

      end #if
    end #j
  end #i

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

  #Connectance before excluding 'made' things
  Sorig = length(find(x->x=='n',diag(int_m)));
  Lorig = sum(t_mpre)/2;
  Corig = Lorig/(Sorig^2);

  #2: if player A contains an "m" with player B, player B is "i" with everything
  # except if it is "n" with player A
  #Which species 'make things?'
  #NOTE: multiple A,B,C can make the same D

  #NOTE: Should 'made' things only be needed and not assimilated, which implies biomass flow? i.e. trophic int?
  madev = Array{Int64}(0);
  for i = 2:num_play
    int_v = copy(int_m[i,:]);
    made = find(x->x=='m',int_v);
    l_made = length(made)
    if l_made > 0
      for k = 1:l_made
        #Compile list of made things
        push!(madev,made[k])
        #The made thing ignores everything... including diag
        int_m[made[k],:] = 'i'
        #Update trophic matrix:
        #the made thing eats nothing
        #It is not a living thing, so nothing 'eats' it
        t_mpre[made[k],:] = 0;
        t_mpre[:,made[k]] = 0;
        #what thing(s) make it?
        makers = find(x->x=='m',int_m[:,made[k]])
        int_m[made[k],makers] = 'n' #Except the thing that makes it
      end
    end
  end

  # #MANY I-I rows/columns? nope
  # itestrc = Array{Bool}(num_play);
  # itestr = Array{Bool}(num_play);
  # itestc = Array{Bool}(num_play);
  # for i=1:num_play
  #   itestr[i] = all(int_m[i,:].=='i');
  #   itestc[i] = all(int_m[:,i].=='i');
  #   itestrc[i] = (itestr[i]== true && itestc[i]==true);
  # end

  #Which rows/columns are not animate objects?
  sprc = find(x->x=='n',diag(int_m));
  num_sp=length(sprc);
  sp_m = Array{Char}(num_sp+1,num_sp+1);
  sp_m[1,:] = hcat('i',int_m[1,sprc]);
  sp_m[:,1] = vcat('i',int_m[sprc,1]);
  sp_m[collect(2:num_sp+1),collect(2:num_sp+1)] = int_m[sprc,sprc];

  #How many 'made things?'
  Snew = length(find(x->x=='n',diag(int_m)));
  Lnew = sum(t_mpre)/2;
  Cnew = Lnew/(Snew^2);

  #Derive the trophic adjacency matrix, which includes all 'ai','aa','an' interactions
  # t_m = zeros(num_play,num_play);
  # for i=1:num_play
  #   for j=1:num_play
  #     if i > j
  #       if int_m[i,j] == 'a' && int_m[j,i] == 'i'
  #         t_m[i,j] = 1
  #         t_m[j,i] = 1
  #       end
  #       if int_m[i,j] == 'i' && int_m[j,i] == 'a'
  #         t_m[i,j] = 1
  #         t_m[j,i] = 1
  #       end
  #       if int_m[i,j] == 'a' && int_m[j,i] == 'n'
  #         t_m[i,j] = 1
  #         t_m[j,i] = 1
  #       end
  #       if int_m[i,j] == 'n' && int_m[j,i] == 'a'
  #         t_m[i,j] = 1
  #         t_m[j,i] = 1
  #       end
  #       if int_m[i,j] == 'a' && int_m[j,i] == 'a'
  #         t_m[i,j] = 1
  #         t_m[j,i] = 1
  #       end
  #     end
  #   end
  # end


return(int_m, sp_m, t_mpre)

end
