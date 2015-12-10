function build_template_degrees(num_play, pw_prob, tr_avoid)

  num_play = 50

  p_n=0.02
  p_e=0.05
  p_m=0.1
  p_a=0
  p_i= 1 - sum([p_n,p_e,p_m,p_a]) #Ignore with 1 - pr(sum(other))

  #Defining paiwise probabilities
  pw_prob = [
    pr_ne = p_n*(p_e/(p_e+p_n+p_i+p_m)),
    pr_nn = p_n*(p_n/(p_e+p_n+p_i+p_m)),
    pr_ni = p_n*(p_i/(p_e+p_n+p_i+p_m)),
    pr_nm = p_n*(p_m/(p_e+p_n+p_i+p_m)),
    pr_ia = p_i*(p_a/(p_e+p_a+p_n+p_i)),
    pr_ie = p_i*(p_e/(p_e+p_a+p_n+p_i)),
    pr_ii = p_i*(p_i/(p_e+p_a+p_n+p_i)),
    pr_aa = p_a*(p_a/(p_a+p_i)),
    pr_ee = p_e*(p_e/(p_i+p_n+p_e))
  ]



  #Fill out matrix based on probabilities above
  prob_line = cumsum(sort(pw_prob))

  # Distributing trophic interactions (e's) according to a degree distribution
  #probabilities are for the whole matrix
  #this means we should have N.e<-p.e*(num_play*(num_play-1)) trophic interactions (e)
  #thus the mean degree should be mean.k=N.e/num_play
  #we can sample degrre from a exponential ditribution with mean (1/rate) equal to mean.k; degrees=rexp(num_play,1/mean.k)


  #Still, the degree distrobution scales with fw size.
  N_e = p_e*(num_play*(num_play-1))

  mean_k = N_e/(num_play-1)

  aux=1
  expdist = Exponential(1/mean_k)
  while aux==1
    degrees = rand(expdist,num_play-1)
    degrees = round(Int64,degrees)
    degrees = [0;degrees] #degree of the sun is 0
    if maximum(degrees)<num_play-1
      aux = 0
    end
  end

  #Create an empty character array with dimensions equal to the number of players
  int_m = Array(Char,num_play,num_play)
  #Set array equal to zero
  int_m[1:num_play*num_play] = '0'

  for i = 2:num_play

    #Assigning initial trophic interactions
    vec = [1:num_play]
    #Remove ith element
    deleteat!(vec,i)
    #What is the degree of player i?
    k = degrees[i]
    resource = rand(vec,k)
    int_m[i,[resource]] = 'e'

    #Assigning complimentary interactions
    #First disconsider ee interactions
    ee_int = find(x -> x == 'e', int_m[[resource],i])

    if length(ee_int) >= 1
      #Remove elements ee_int
      deleteat!(resource,ee_int)
    end

    #randomly determine ei or en (ee already emerge from assigning e's)
    bindist = Binomial(1,pw_prob[1]/pw_prob[6]) #PROBLEMS NEVER DRAWS ONES
    bin_draw = rand(bindist,length(resource))

    #Faculatative mutualism
    int_m[find(x->x==1,resource.*bin_draw),i]<-"n" #assigning n's according to random draw
    #Asymmetric predation
    int_m[find(x->x==0,resource.*bin_draw),i]<-"i"#assigning i's according to random draw

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
    for k = 1:l_made
