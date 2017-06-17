# L = zeros(6,6);
# L[1:36] = [0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
# L = transpose(L);


function trophicalc(L)
  s = size(L)[1];
  #Find primary producers
  rowsums = [sum(L[i,:]) for i=1:s];
  leaves = find(x->x==0,rowsums);
  branches = collect(1:s);
  deleteat!(branches,leaves);
  
  maxHeight = s - length(leaves) - 1;
  index = collect(1.0:maxHeight);
  paths = Array(Array{Float64},length(index));
  pathSteps = Array(Array{Float64},length(index));
  for i=1:length(index)
    paths[i] = L^index[i];
    pathSteps[i] = (index[i])*paths[i];
  end
  paths=reduce(+,paths)[branches,leaves];
  pathSteps=reduce(+,pathSteps)[branches,leaves];
  rowsumspaths = zeros(Float64,s);
  rowsumspathSteps=zeros(Float64,s);
  #Calc rowsums over internal species
  for i=1:length(branches)
    rowsumspaths[i] = sum(paths[i,:]);
    rowsumspathSteps[i] = sum(pathSteps[i,:]);
  end
  branchtrophic = (rowsumspathSteps./rowsumspaths);
  out = zeros(Float64,s);
  tic = 0;
  for i=branches
    tic = tic+1;
    out[i]=branchtrophic[tic];
  end
  return(out);
end
