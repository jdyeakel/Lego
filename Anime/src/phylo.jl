function phylo(int_m)
    
    ml = size(int_m)[1];
    
    #Compute the distance metric
    dm = Array{Float64}(ml,ml);
    
    for i=1:ml
        for j=1:ml
            
            simvec = int_m[i,:] .== int_m[j,:];
            
            dist = 1 - (sum(simvec)/length(simvec))
            
            dm[i,j] = dist
        end
    end
    
    return dm
    
end
    