function trophicwidth(L)
  sl = size(L)[1];
  rowsums = Array{Float64}(sl);
  [rowsums[i]=sum(L[i,:]) for i=1:sl];
  keep = [1;find(x->x>0,rowsums)];
  culledmatrix = L[keep,keep];
  
  lm = length(keep);
  
  colsums = Array{Float64}(lm);
  [colsums[i]=sum(culledmatrix[:,i]) for i=1:lm];
  
  
  tw = zeros(Float64,lm);
  # consload = zeros(Float64,lm);
  for i=2:lm
    #find resources for consumer i
    res = find(x->x==1,culledmatrix[i,2:lm]); #don't count sharing basal resource
    #The proportion of consumers that eat consumer i's resources
    consload = colsums[res]/(lm-1); #lm is number of species, -1 to discount basal resource
    #mean proportion of consumers that eat consumer i's resources
    tw[i] = mean(consload);
  end
  #get rid of nan
  tw[isnan(tw)]=0;
  
  return(tw);

end
