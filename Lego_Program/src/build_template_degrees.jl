function build_template_degrees(num_play, probs)


  #Defining paiwise probabilities
  #These are probabilities of pairwise interactions within the whole universe of possible interactions (15)
  #THIS IS PROBABLY WRONG
  # pw_prob_init = [
  #   pr_na = p_n*(p_a/(p_a+p_n+p_i+p_m)),
  #   pr_nn = p_n*(p_n/(p_a+p_n+p_i+p_m)),
  #   pr_ni = p_n*(p_i/(p_a+p_n+p_i+p_m)),
  #   pr_nm = p_n*(p_m/(p_a+p_n+p_i+p_m)),
  #   pr_ie = p_i*(p_e/(p_a+p_e+p_n+p_i)),
  #   pr_ia = p_i*(p_a/(p_a+p_e+p_n+p_i)),
  #   pr_ii = p_i*(p_i/(p_a+p_e+p_n+p_i)),
  #   pr_ee = p_e*(p_e/(p_e+p_i)),
  #   pr_aa = p_a*(p_a/(p_i+p_n+p_a))
  # ]

  #THIS IS PROBABLY RIGHT
  pw_prob_init = [
    pr_na = p_n*(p_a/(p_a+p_n+p_i+p_m)) + p_a*(p_n/(p_a+p_i+p_n)),
    pr_nn = p_n*(p_n/(p_a+p_n+p_i+p_m)),
    pr_ni = p_n*(p_i/(p_a+p_n+p_i+p_m)) + p_i*(p_n/(p_a+p_n+p_i+p_e)),
    pr_nm = p_n*(p_m/(p_a+p_n+p_i+p_m)) + p_m*(p_n/p_n),
    pr_ie = p_i*(p_e/(p_a+p_e+p_n+p_i)) + p_e*(p_i/(p_i+p_e)),
    pr_ia = p_i*(p_a/(p_a+p_e+p_n+p_i)) + p_a*(p_i/(p_a+p_i+p_n)),
    pr_ii = p_i*(p_i/(p_a+p_e+p_n+p_i)),
    pr_ee = p_e*(p_e/(p_e+p_i)),
    pr_aa = p_a*(p_a/(p_i+p_n+p_a))
  ];


  # pw_prob = copy(pw_prob_init) / sum(pw_prob_init)

  #Fill out matrix based on probabilities above
  prob_line = cumsum(sort(pw_prob_init));

  # Distributing trophic interactions (e's) according to a degree distribution
  #probabilities are for the whole matrix
  #this means we should have N.e<-p.e*(num_play*(num_play-1)) trophic interactions (e)
  #thus the mean degree should be mean.k=N.e/num_play
  #we can sample degrre from a exponential ditribution with mean (1/rate) equal to mean.k; degrees=rexp(num_play,1/mean.k)


  #Still, the degree distrobution scales with fw size.
  #Expected number of trophic links
  N_a = p_a*(num_play*(num_play-1));
  #Expected number of trophic links per species (discounting cannibalism)
  mean_k = N_a/(num_play-1);

  #Assuming trophic links follow an exponential distribution
  #Draw number of trophic links per node, with element 1=sun (k=0)
  expdist = Exponential(1/mean_k);
  degrees = Array(Int64);
  aux=1
  while aux == 1
    degrees = rand(expdist,num_play-1);
    degrees = round(Int64,degrees);
    degrees = [0;copy(degrees)]; #degree of the sun is 0
    if maximum(degrees)<num_play-1
      aux = 0;
    end
  end

  #Create an empty character array with dimensions equal to the number of players
  int_m = Array(Char,num_play,num_play)
  #Set array equal to zero
  int_m[1:num_play*num_play] = '0';

  for i = 2:num_play

    #Assigning initial trophic interactions
    vec = collect(1:num_play);
    #Remove ith element
    deleteat!(vec,i)
    #What is the degree of player i?
    k = degrees[i]
    #Randomly choose the identities of prey
    resource = rand(vec,k)
    #Establish these prey in the interaction matrix
    int_m[i,collect(resource)] = 'a'
    #Note: to vectorize a row from int_m (so that it is an Array{Char,1}), we would write v = int_m[i,:][:]

    #Assigning complimentary interactions
    #First disconsider aa interactions
    aa_int = find(x -> x == 'a', int_m[collect(resource),i])

    if length(aa_int) >= 1
      #Remove elements ee_int
      deleteat!(resource,aa_int)
    end

    #randomly determine ai or an (aa already emerges from assigning a's)
    #Draw from a binomial: which interaciton will be Asymmetric predation?
    #Draw 1 = Asymmetric predation; Draw 0 = facultative mutualism
    pr_i_given_a = (p_i/(p_i+p_n));
    bindist = Binomial(1,pr_i_given_a) #PROBLEMS NEVER DRAWS ONES
    bin_draw = rand(bindist,length(resource))

    #Faculatative mutualism
    res_i = resource[find(x->x==1,bin_draw)]
    int_m[[res_i],i] = 'i' #assigning n's according to random draw
    #Asymmetric predation
    res_n = resource[find(x->x==0,bin_draw)]
    int_m[[res_n],i] = 'n' #assigning i's according to random draw

  end

  int_labels = ["ne","nn","ni","nm","ia","ie","ii","aa","ee"]

  pw_prob_new = pw_prob
  deleteat!(pw_prob_new,[1,6,9])
  new_int_labels = int_labels
  pw_prob_new = pw_prob_new/sum(pw_prob_new)
  deleteat!(new_int_labels,[1,6,9])

  prob_line = cumsum(sort(pw_prob_new))



  for i = 2:num_play
    for j = 2:num_play
      #Only choose interactions for empty elements (int_m = '0')
      #Only do this for the triangular part of the matrix
      if int_m[i,j] == '0' && i < j

        r_draw = rand()

        #N:N
        if r_draw < prob_line[3]
          index = find(x -> x == "nn",new_int_labels)
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
        if r_draw > prob_line[3] && r_draw < prob_line[4]
          index = find(x -> x == "nm",new_int_labels)
          mij = Array(Char,2)
          mij[1] = rand([new_int_labels[index][1],new_int_labels[index][2]])
          if mij[1] == new_int_labels[index][1]
            mij[2] = new_int_labels[index][2]
          else mij[2] = new_int_labels[index][1]
          end
          int_m[i,j] = mij[1]
          int_m[j,i] = mij[2]
        end

        #N:I
        if r_draw > prob_line[4] && r_draw < prob_line[5]
          index = find(x -> x == "ni",new_int_labels)
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
        if r_draw > prob_line[5] && r_draw < prob_line[6]
          index = find(x -> x == "ii",new_int_labels)
          mij = Array(Char,2)
          mij[1] = rand([new_int_labels[index][1][1],new_int_labels[index][1][2]])
          if mij[1] == new_int_labels[index][1][1]
            mij[2] = new_int_labels[index][1][2]
          else mij[2] = new_int_labels[index][1][1]
          end
          int_m[i,j] = mij[1]
          int_m[j,i] = mij[2]
        end

      end #if
    end #j
  end #i

  #Filling in the diagonal
  #Determine which are active players and which aren't
  for i = 1:num_play
    if length(find(x->x=='e',int_m[i,:])) > 0
      int_m[i,i] = 'n'
    else int_m[i,i] = 'i'
    end
  end


  #Implement posthoc Rules

  #1: Row/Col 1 is the sun
  int_m[1,:] = 'i'
  #All column elements that ARENT 'e' are 'ignore'
  int_m[find(x->x!='e',int_m[:,1])] = 'i'

  #2: if player A contains an "m" with player B, player B is "i" with everything
  # except if it is "n" with player A
  #Which species 'make things?'
  #NOTE: multiple A,B,C can make the same D

  for i = 1:num_play
    int_v = int_m[i,:]
    made = find(x->x=='m',int_v)
    l_made = length(made)
    if l_made > 0
      for k = 1:l_made
        #The made thing ignores everything...
        int_m[made[k],:] = 'i'
        #what thing(s) make it?
        makers = find(x->x=='m',int_m[:,made[k]])
        int_m[made[k],makers] = 'n' #Except the thing that makes it
      end
    end
  end

return(int_m)

end
