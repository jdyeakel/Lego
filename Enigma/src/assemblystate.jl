function assemblystate(S,probs,lambda)


    MaxN = convert(Int64,floor(S + S*lambda));

    int_m, tp_m, tind_m, mp_m, mind_m = intmatrixv3(S,lambda,probs);

    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id = preamble_defs(int_m);

    N = size(int_m)[1];

    #Define all possible states
    states = collect(combinations(collect(1:N)));
    states = states[N+1:length(states)];
    states = states[in.(1,states)];
    #how many states
    nstates = length(states);

    #Determine impossible states
    passtest = zeros(Int64,nstates) .+ 1;
    keepstates = Array{Int64}(undef,0);
    for i=1:nstates
        speciesobjects = states[i];
        adjacencymatrix = a_b[speciesobjects,speciesobjects] .+ m_b[speciesobjects,speciesobjects]';
        g = DiGraph(adjacencymatrix');
        paths = gdistances(g,1)
        if maximum(paths) < N+1;
            passtest[i] *= 1;
        else
            passtest[i] *= 0;
        end
        
        #NOTE: take out states without complete set of species/object pairs
        observedspecies = speciesobjects[speciesobjects .<= S];
        expectedobjects = findall(!iszero,vec(sum(m_b[observedspecies,:],dims=1)));
        
        observedobjects = setdiff(speciesobjects,observedspecies);
        expectedspecies = findall(!iszero,vec(sum(m_b[:,observedobjects],dims=2)));
        
        if observedspecies == expectedspecies || observedobjects == expectedobjects
            passtest[i] *= 1;
        else
            passtest[i] *= 0;
        end
        
        
    end

    possiblestates = states[findall(!iszero,passtest)];
    statespace = sum(passtest)/nstates
    lstate = length(possiblestates);

    # SparseArray
    transm = spzeros(lstate,lstate);
    for i=1:lstate
        # print(string(i,'_'))
        statei = copy(possiblestates[i]);
        deleteat!(statei,1)
        colonizers = potcol(sp_v,int_id,statei,a_b,n_b0,0,1);
        # newstates = Array{Array}(undef,length(colonizers));
        newstatesloc = Array{Int64}(undef,length(colonizers));
        for j=1:length(colonizers)
            newstates = sort([1;statei;colonizers[j]]);
            newstatesloc[j] = findall(x->x==newstates,possiblestates)[1];
        end
        transm[i,newstatesloc] .= 1.0;
    end

    primarystates = findall(x->x==2,length.(possiblestates));
    transg = DiGraph(transm);
    trimlist = zeros(Int64,lstate,length(primarystates));
    for i=1:length(primarystates)
        paths = gdistances(transg,primarystates[i]);
        tokeep = findall(x->x<(S+1),paths);
        trimlist[tokeep,i] .= 1;
    end 
    connectedstates = findall(!iszero,vec(sum(trimlist,dims=2)));
    transmconnected = transm[connectedstates,connectedstates];
    lstateconnected = length(connectedstates);
    
    return(int_m,transmconnected,connectedstates,possiblestates) 
    
end
