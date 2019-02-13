function nichemodelweb(S,C) 
    
    
    #list of random niche values
    n = sort(rand(S));
    
    #Define beta distribution
    a = 1;
    b = (1/(2*C)) - 1;
    bdist = Beta(a,b);
    
    #range
    r = n.*rand(bdist,S);
    
    #center of range
    c = rand.(Uniform.(r/2,n));
    
    #Find species that fall within the range of each other
    prey = Array{Array{Int64}}(undef,S);
    [prey[i] = findall(x-> (x > c[i] - (r[i]/2) && x < c[i] + (r[i]/2)), n) for i=1:S];
    

    adjmatrix = zeros(Bool,S,S);
		for i=1:S
    	adjmatrix[i,prey[i]] .= true;
		end
			
    
    ladj = size(adjmatrix)[1];
    keep = findall(!iszero,vec(sum(adjmatrix,dims=2))+vec(sum(adjmatrix,dims=1)));
    niche = n;
    
    
    while length(keep) < ladj
        ladj = size(adjmatrix)[1];
        keep = findall(!iszero,vec(sum(adjmatrix,dims=2))+vec(sum(adjmatrix,dims=1)));
        adjmatrix = adjmatrix[keep,keep];
        niche = niche[keep];
    end
    
    
    
    return adjmatrix, niche
    
end

# R"""
# image($(adjmatrix),col=grey(c(0,1)))
# """
