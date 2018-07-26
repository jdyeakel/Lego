function assembly2(int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
    athresh,nthresh,tmax)
    
    N = size(int_m)[1];
    v = collect(1:N);
    cid = falses(N);
    
    
    S = length(sp_v) + 1;
    vsp = collect(1:S);
    spcid = falses(S);
    
    vob = collect(S+1:N);
    O = length(vob);
    ocid = falses(O);
    
    sprich = Array{Int64}(0);
    rich = Array{Int64}(0;)
    clock = Array{Float64}(0);
    
    #Build the strength matrix apriori
    strength = (pi*sum(n_b0,2)) .- sum(a_b,2);
    smatrix = Array{Float64}(copy(a_b));
    for i=1:size(a_b)[1]
        smatrix[i,:] = a_b[i,:] * strength[i];
    end
    # smatrix[find(iszero,a_b)] = NaN;
    
    
    t=0;
    it = 0;
    while t < tmax
        
        cid_old = copy(cid);
        #Which are species?
        # spcid = intersect(sp_v,v[cid]);
        #Which are objects?
        # ocid = setdiff(v[cid],spcid);
        
        
        #COUNT POTENTIAL COLONIZERS
        trophiclinked = setdiff(int_id[(sum(a_b[:,[1;v[cid]]],2) .> 0)[:,1]],v[cid]);
        #For each trophiclinked, count number of assimilate and need interactions in system
        #Determine in the proportion that already exists is >= the threshold
        a_fill = ((sum(a_b[trophiclinked,[1;v[cid]]],2)./sum(a_b[trophiclinked,:],2)) .>= athresh)[:,1];
        prop_n = sum(n_b0[trophiclinked,[1;v[cid]]],2)./sum(n_b0[trophiclinked,:],2);
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
        a_fill = ((sum(a_b[vsp[spcid],[1;v[cid]]],2)./sum(a_b[vsp[spcid],:],2)) .> athresh)[:,1];
        prop_n = sum(n_b0[vsp[spcid],[1;v[cid]]],2)./sum(n_b0[vsp[spcid],:],2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] = 1;
        n_fill = (prop_n .>= nthresh)[:,1];
        #Build a list of those that pass each individual test
        a_pass = vsp[spcid][a_fill];
        n_pass = vsp[spcid][n_fill];
        
        #which species pass?
        survivors = intersect(a_pass,n_pass);
        # spext1 = setdiff(vsp[spcid],survivors);
        
        vspext = falses(S);
        vspext[vsp[spcid]] = true;
        vspext[survivors] = false;
        #Alternative approach
        
        #define subset of smatrix for community at this timestep
        cmatrix = smatrix[vsp[spcid],v[cid]];
        #define max values for community at this timestep
        cmax = vec(findmax(cmatrix,1)[1]);
        
        # NOTE: THis is not faster
        # cinv = 1./cmax;
        # m = (strength[spcid] * cinv');
        # m0 = m .* a_b[spcid,cid];
        # spext2 = spcid[vec(findmax(m0,2)[1] .< 1)];
        # 
        prext_comp = trues(length(vsp[spcid]));
        for i=1:length(vsp[spcid])
            ieats = (a_b[i,v[cid]] .== true);
            if any(ieats)
                if any(strength[vsp[spcid][i]] .>= cmax[ieats]); #cmax[find(!isnan,cmatrix[i,:])]);
                    #species is NOT added to pool
                    prext_comp[i] = false;
                end
            else
                prext_comp[i] = false;
            end
        end
        spext2 = vsp[spcid][prext_comp];
        # 
        # spext2 = Array{Int64}(0);
        # spext = unique([spext1;spext2]);
        vspext[spext2] = true;
        lspext = sum(vspext);
        spext = vsp[vspext];
        
        
        #COUNT POTENTIAL EXTINCT OBJECTS
        obext = vob[ocid][find(iszero,sum(m_b[vsp[spcid],vob[ocid]],1))];
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
            c_ob = setdiff(c_ob_full,v[cid_old]);
            colonize = [c_sp;c_ob];
            #Update community
            cid[colonize] = true;
            spcid[c_sp] = true;
            ocid[c_ob-S] = true;
        end
        
        if re > (lcol/levents) && re < ((lcol + lspext)/levents)
            
            #SPECIES EXTINCTION FUNCTION
            #select species to go extinct
            sp_bye = rand(spext,1);
            cid[sp_bye] = false;
            spcid[sp_bye] = false;
        end
        
        if re > ((lcol + lspext)/levents)
            
            #OBJECT EXTINCTION FUNCTION
            ob_bye = rand(obext,1);
            cid[ob_bye] = false;
            ocid[ob_bye-S] = false;
            
        end
        
        #NOTE - updating CID....
        #CID[t,cid] = true;
        push!(sprich,sum(spcid));
        push!(rich,sum(cid));
    end #end time steps
    
    return(
    sprich,
    rich,
    clock
    )
    
end
