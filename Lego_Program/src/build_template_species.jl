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
  num_play = Int64(floor(S + (S^2)*pr_nm));

  #Generate interaction template
  int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_degrees(num_play,probs,ppweight);

  S_real = length(find(x->x=='n',diag(int_m)));
  
  error = 5;
  S_acceptable = collect(400-error:400+error);
  
  howlong = 1;
  while in(S_real,S_acceptable) == false
    howlong = howlong + 1;
    S_diff = S - S_real;
    plusminus = sign(S_diff);
    dhop = Poisson(abs(S_diff));
    np_hop = rand(dhop);

    num_play = num_play+(np_hop*plusminus);

    #Generate interaction template
    int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_degrees(num_play,probs,ppweight);

    S_real = length(find(x->x=='n',diag(int_m)));
  end
  # println("how long? ",howlong)  

return(int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue)

end
