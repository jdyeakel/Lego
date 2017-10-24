function Jmatrix(m_tp_full,m_mp_full)
    
    
    #In this function, we take 1) recorded trophic interactions, and 2) recorded mutualistic interactions, and assign (+1,-1) and (+1,+1) to the Jacobian, respectively
    
    #eliminate non-existant species (for un-full communities)
    rowsums = zeros(size(m_tp_full,1));
    [rowsums[i] = sum(m_tp_full[i,:]) for i=1:size(m_tp_full,1)];
    keep = find(x->x>0,rowsums);
    m_tp = m_tp_full[keep,keep];
    m_mp = m_mp_full[keep,keep];
    
    #matrix size
    lm = size(m_tp,1);
    
    #create empty interaction matrix
    Jmatrix = zeros(Int64,lm,lm);
    #set diagonal
    for i=1:lm
        Jmatrix[i,i] = -1;
    end
    
    #find predation (+,-) interactions
    tint = find(x->x==1,m_tp);
    tintloc = [ind2sub(size(m_tp), i) for i in tint];
    #find mutualistic (+,+) interactions
    mint = find(x->x==1,m_mp);
    mintloc = [ind2sub(size(m_mp), i) for i in mint];
    
    #Assign predator-prey interactions (+1,-1)
    for i=1:length(tintloc)
        Jmatrix[tintloc[i][1],tintloc[i][2]] = 1;
        Jmatrix[tintloc[i][2],tintloc[i][1]] = -1;
    end
    #Assign mutualistic interactions (+1,+1)
    for i=1:length(mintloc)
        Jmatrix[mintloc[i][1],mintloc[i][2]] = 1;
        Jmatrix[mintloc[i][2],mintloc[i][1]] = 1;
    end
    
    #trim off non-interactions
    #which species have no interactions?
    # colsums = Array{Float64}(lm);
    # [colsums[i]=sum(culledmatrix[:,i]) for i=1:lm];
    # find(x->x==0,sum())
    # 
    
    
    
    
    return(Jmatrix)
    
    
end

