function extinct_func(cid,c_m,crev_m,com_sparse,com_tp,com_tind)

  #Prior community structure
  cidold = copy(cid);
  c_mold = copy(c_m);
  crev_mold = copy(crev_m);


  #Randomly select the species to eliminate
  esp = rand(cid);

  #Determine the cascade of consequent extinctions

  #1) What does it make? eliminate those things
  em = find(x->x=='m',int_m[esp,:]);


end
