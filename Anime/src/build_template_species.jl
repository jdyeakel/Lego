function build_template_species(S, probs, ppweight)

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


  #Initial guess
  # num_play = Int64(floor((S^2)*pr_nm));
  
  num_play = Int64(floor(S + (S^2)*pr_nm));
  np = Array{Int64}(0);
  push!(np,num_play);
  #num_play = convert(Int64,round(num_play/2,0));
  
  # @label p;
  
  #Generate interaction template
  int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = intmatrix(num_play,probs,ppweight);

  S_real = length(find(x->x=='n',diag(int_m)));

  error = 10;
  S_acceptable = collect(S-error:S+error);

  howlong = 1;
  sdiffvec = Array{Int64}(0);
  
  et = 0;
  tic();
  while in(S_real,S_acceptable) == false
    # #This breaks the while loop and creates a new int_m 
    # et = toc();
    # if et > 10
    #     @goto p
    # end
    howlong = howlong + 1;
    S_diff = S - S_real;
    push!(sdiffvec,S_diff);
    plusminus = sign(S_diff);
    
    #Choose the next step size for the size of int_m
    dhop = Poisson(2*abs(S_diff));
    # np_hop = rand(dhop);
    np_hop = convert(Int64,round(rand(dhop)*(S/S_real),0));

    num_play = num_play+(np_hop*plusminus);
    push!(np,num_play);

    #Generate interaction template
    int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m = intmatrix(num_play,probs,ppweight);

    S_real = length(find(x->x=='n',diag(int_m)));
    println("S_diff = ",S_diff)
  end
  # println("how long? ",howlong)

return(int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m)

end
