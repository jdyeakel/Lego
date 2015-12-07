function rexp(n,mean)
  if mean <= 0.0
        error("mean must be positive")
  end
  out = zeros(n)
  for i = 1:n
    out[i] = -mean*log(rand())
  end
  return(out)
end
