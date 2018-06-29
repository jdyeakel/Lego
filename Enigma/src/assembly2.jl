function assembly2(int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
    athresh,nthresh,tmax)
    
    S = length(sp_v) + 1;
    N = size(int_m)[1];
    # MaxN = convert(Int64,floor(S + S*lambda));
    cid = Array{Int64}(0);
    sprich = Array{Int64}(0);
    rich = Array{Int64}(0;)
    clock = Array{Float64}(0);
    
    #Build the strength matrix apriori
    strength = (pi*sum(n_b0,2)) .- sum(a_b,2);
    smatrix = Array{Float64}(copy(a_b));
    for i=1:size(a_b)[1]
        smatrix[i,:] = a_b[i,:] * strength[i];
    end
    smatrix[find(iszero,a_b)] = NaN;
    
    spvec = collect(2:S);
    obvec = collect(S+1:N);
    spbk = 0;
    obbk = 0;
    t=0;
    it = 0;
    @time while t < tmax
        
        cid_old = copy(cid);
        #Which are species?
        spcid = spvec[1:spbk]; #surprisingly this works with spbk=0
        # spcid = intersect(sp_v,cid);
        #Which are objects?
        ocid = obvec[1:obbk];
        # ocid = setdiff(cid,spcid);
        cid = [spcid;ocid];
        
        #COUNT POTENTIAL COLONIZERS
        nspcid = spvec[spbk+1:N];
        trophiclinked = int_id[vec(sum(a_b[:,[1;cid]],2) .> 0)];
        # trophiclinked = setdiff(int_id[(sum(a_b[:,[1;cid]],2) .> 0)[:,1]],cid);
        #For each trophiclinked, count number of assimilate and need interactions in system
        #Determine in the proportion that already exists is >= the threshold
        a_fill = ((sum(a_b[trophiclinked,[1;cid]],2)./sum(a_b[trophiclinked,:],2)) .>= athresh)[:,1];
        prop_n = sum(n_b0[trophiclinked,[1;cid]],2)./sum(n_b0[trophiclinked,:],2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] = 1;
        n_fill = (prop_n .>= nthresh)[:,1];
        #Build a list of those that pass each individual test
        a_pass = trophiclinked[a_fill];
        n_pass = trophiclinked[n_fill];
        #Build a list of those that pass both tests (potential colonizers)
        col = intersect(a_pass,n_pass);
        #Count the number that pass
        lcol = length(col);
        
        
        #COUNT POTENTIAL EXTINCT SPECIES
        #1) By not fulfilling Eat/Need thresholds
        #Re-calculate athresh and nthresh (> athresh; >= nthresh)
        a_fill = ((sum(a_b[spcid,[1;cid]],2)./sum(a_b[spcid,:],2)) .> athresh)[:,1];
        prop_n = sum(n_b0[spcid,[1;cid]],2)./sum(n_b0[spcid,:],2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] = 1;
        n_fill = (prop_n .>= nthresh)[:,1];
        #Build a list of those that pass each individual test
        a_pass = spcid[a_fill];
        n_pass = spcid[n_fill];
        #which species pass?
        survivors = intersect(a_pass,n_pass);
        spext1 = setdiff(spcid,survivors);
        
        #2) By not having the max strength for consuming at least one resource
        # strength = sum(n_b0[spcid,[1;cid]],2) .- sum(a_b[spcid,[1;cid]],2);
        # smatrix = Array{Float64}(copy(a_b[spcid,[1;cid]]));
        # for i=1:length(spcid)
        #     smatrix[i,:] = a_b[spcid[i],[1;cid]] * strength[i];
        # end
        # smatrix[find(iszero,a_b[spcid,[1;cid]])] = NaN;
        
        # #NOTE: This is probably the slowest part
        # prext_comp = Array{Bool}(length(spcid));
        # for i=1:length(spcid)
        #     cmatrix = smatrix[:,find(!iszero,a_b[spcid[i],[1;cid]])];
        #     #Proportion of strength-max foods
        #     propmax = sum(cmatrix[i,:] .>= findmax(cmatrix,1)[1]')/length(cmatrix[i,:]);
        #     if propmax > 0
        #         prext_comp[i] = false;
        #     else
        #         prext_comp[i] = true;
        #     end
        # end
        # spext2 = spcid[prext_comp];
        
        #Alternative approach
        prext_comp = Array{Bool}(length(spcid));
        #define subset of smatrix for community at this timestep
        cmatrix = smatrix[spcid,cid];
        #define max values for community at this timestep
        cmax = vec(findmax(cmatrix,1)[1]);
        for i=1:length(spcid)
            if any(strength[spcid[i]] .>= cmax[find(!isnan,cmatrix[i,:])]);
                #species is NOT added to pool
                prext_comp[i] = false;
            else 
                #species IS added to pool
                prext_comp[i] = true;
            end
        end
        spext2 = spcid[prext_comp];
        
        spext = unique([spext1;spext2]);
        lspext = length(spext);
        
        
        #COUNT POTENTIAL EXTINCT OBJECTS
        obext = ocid[find(iszero,sum(m_b[spcid,ocid],1))];
        lobext = length(obext);
        
        levents = sum([lcol;lspext;lobext]);
        
        dt = 1/levents;
        
        t += dt;
        it += 1;
        push!(clock,t);
        
        #Choose a random event
        re = rand();
        
        if re < (lcol/levents)
            
            #COLONIZATION FUNCTION
            c_sp = rand(col);
            #select made objects as well
            c_ob_full = find(isodd,m_b[c_sp,:]);
            #Only bring in objects that aren't already in the community
            c_ob = setdiff(c_ob_full,cid_old);
            colonize = [c_sp;c_ob];
            #Update community
            cid = [cid_old;colonize];
            
        end
        
        if re > (lcol/levents) && re < ((lcol + lspext)/levents)
            
            #SPECIES EXTINCTION FUNCTION
            #select species to go extinct
            sp_bye = rand(spext,1);
            cid = setdiff(cid_old,sp_bye);
            
        end
        
        if re > ((lcol + lspext)/levents)
            
            #OBJECT EXTINCTION FUNCTION
            ob_bye = rand(obext,1);
            cid = setdiff(cid_old,ob_bye);
            
        end
        
        #NOTE - updating CID....
        #CID[t,cid] = true;
        push!(sprich,length(spcid));
        push!(rich,length(cid));
    end #end time steps
    
    return(
    sprich,
    rich
    )
    
end
