function build_template_degrees(num_players, pw_prob, tr_avoid)

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
  while aux==1
    expdist = Exponential(1/mean_k)
    degrees = rexp(expdist,num_play-1)
    degrees = rand(Int64,degrees)
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

    if size(ee_int)[1] >= 1
      #Remove elements ee_int
      deleteat!(resource,ee_int)

      #randomly determine ei or en (ee already emerge from assigning e's)
      bindist = Binomial(1,pw_prob[1]/pw_prob[6])
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
        if int_m[i,j] == '0'
          if i < j

            r_draw = rand()

            #N:N
            if r_draw < prob_line[3]
              mij = Array(Char,2)
              mij[1] = rand([new_int_labels[1][1],new_int_labels[1][2]])
              if mij[1] == new_int_labels[1][1]
                mij[2] = new_int_labels[1][2]
              else mij[2] = new_int_labels[1][1]
              end
              int_m[i,j] = mij[1]
              int_m[j,i] = mij[2]
            end

            #N:M
            if r_draw > prob_line[3] && r_draw < prob_line[4]
              mij = Array(Char,2)
              mij[1] = rand([new_int_labels[2][1],new_int_labels[2][2]])
              if mij[1] == new_int_labels[2][1]
                mij[2] = new_int_labels[2][2]
              else mij[2] = new_int_labels[2][1]
              end
              int_m[i,j] = mij[1]
              int_m[j,i] = mij[2]
            end
