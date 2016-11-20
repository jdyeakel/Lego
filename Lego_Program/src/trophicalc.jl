function trophicalc(Slist,cid,com_tp)

  #Identify species and their position within com_tp and/or com_tind (whatever the input)
  lsp = length(cid);
  t_loc = zeros(Int64,lsp);
  for i=1:lsp
    t_loc[i] = find(x->x==cid[i],Slist)[1] + 1;
  end


  



end
