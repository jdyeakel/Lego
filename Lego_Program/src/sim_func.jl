function sim_func(seq1,seq2)

  #determine where both seq1 and seq1 have i,i & delete those interactions
  todelete = Array{Int64}(0);
  for i=1:num_play
    if seq1[i] == 'i' && seq2[i] == 'i'
      push!(todelete,i);
    end
  end
  deleteat!(seq1,todelete);
  deleteat!(seq2,todelete);

  #Calculate similarity between non-i-i interactions
  sim = sum(seq1 .== seq2)/length(seq1);

  return(sim)

end
