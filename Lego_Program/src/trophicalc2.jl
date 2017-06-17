function trophicalc(L)
  rowsums = Array{Float64}(size(L)[1]);
  [rowsums[i]=sum(L[i,:]) for i=1:size(L)[1]];
  keep = [1;find(x->x>0,rowsums)];
  culledmatrix = L[keep,keep];
  R"""
  library(NetIndices)
  library(MASS)
  rtl<-TrophInd(t($culledmatrix))
  """
  @rget rtl;
  tl = Array(rtl[1]);
  return(tl);
end
