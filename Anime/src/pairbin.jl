function pairbin(m1, m2)
    
    m1t = UpperTriangular(m1)' + LowerTriangular(m1);
    m2t = UpperTriangular(m2)' + LowerTriangular(m2);
    
    lt = LowerTriangular(m1t .* m2t);
    
    sym = ((lt + lt') .>= 1)*1;
    
    return(sym)
    
end