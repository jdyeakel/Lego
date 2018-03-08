function pairbin(m1, m2)
    
    #Diagonal must be zero for this function to work
    m1[diagind(m1)]=0;
    m2[diagind(m2)]=0;
    
    if m1 == m2
        #Already-symmetric links (a-a), (i-i), (n-n)
        #Eliminates any links that are not in both m1 and m2
        lt = UpperTriangular(m1)' .* LowerTriangular(m2);
        sym = convert(Array{Bool},(lt + lt'))*1;
    else 
    
        #This combines directed links in upper and lower matrices to make them symmetric
        #Already-symmetric links will be >1
        m1t = UpperTriangular(m1)' + LowerTriangular(m1);
        m2t = UpperTriangular(m2)' + LowerTriangular(m2);
        
        #Eliminates any links that are not in both m1t and m2t
        lt = LowerTriangular(m1t .* m2t);
        
        sym = convert(Array{Bool},(lt + lt'))*1;
    end
    
    return(sym)
    
end