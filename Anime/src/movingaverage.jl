function movingaverage(vec, window)
  l = size(vec);
  avgvec = zeros(Float64,l);
  avgvec = avgvec[1:l[1]-window,:];
  if length(l) > 1
    for i = 1:l[2] #Over different elements of the timeseries
      for t = 1:(l[1]-window)
        avgvec[t,i] = mean(vec[t:t+window,i]);
      end
    end
  else
    for t = 1:(l[1]-window)
      avgvec[t] = mean(vec[t:t+window]);
    end
  end
  
  
  return(avgvec)
end
