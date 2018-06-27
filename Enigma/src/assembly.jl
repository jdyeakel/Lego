function (int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,tmax)
    
    
    
    sprich = Array{Int64}(tmax);
    prim_ext = Array{Int64}(tmax);
    sec_ext = Array{Int64}(tmax);
    status = Array{Int64}(tmax);
    lpot_col = Array{Int64}(tmax);
    CID = Array{Bool}(tmax,MaxN)*false;
    
    S = length(sp_v) + 1;
    MaxN = convert(Int64,floor(S + S*lambda));
    cid = Array{Int64}(0);
    sprich = Array{Int64}(0;)
    clock = Array{Float64}(0);
    t=0;
    it = 0;
    while t < tmax
        
        cid_old = copy(cid);
        #Which are species?
        spcid = intersect(sp_v,cid);
        #Which are objects?
        ocid = setdiff(cid,spcid);
        
        
        #COUNT POTENTIAL COLONIZERS
        trophiclinked = setdiff(int_id[(sum(a_b[:,[1;cid]],2) .> 0)[:,1]],cid);
        #For each trophiclinked, count number of assimilate and need interactions in system
        #Determine in the proportion that already exists is >= the threshold
        a_fill = ((sum(a_b[trophiclinked,[1;cid]],2)./sum(a_b[trophiclinked,:],2)) .>= a_thresh)[:,1];
        prop_n = sum(n_b0[trophiclinked,[1;cid]],2)./sum(n_b0[trophiclinked,:],2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] = 1;
        n_fill = (prop_n .>= n_thresh)[:,1];
        #Build a list of those that pass each individual test
        a_pass = trophiclinked[a_fill];
        n_pass = trophiclinked[n_fill];
        #Build a list of those that pass both tests (potential colonizers)
        col = intersect(a_pass,n_pass);
        #Count the number that pass
        lcol = length(col);
        
        
        #COUNT POTENTIAL EXTINCT SPECIES
        #1) By not fulfilling Eat/Need thresholds
        #Re-calculate a_thresh and n_thresh (> a_thresh; >= n_thresh)
        a_fill = ((sum(a_b[spcid,[1;cid]],2)./sum(a_b[spcid,:],2)) .> a_thresh)[:,1];
        prop_n = sum(n_b0[spcid,[1;cid]],2)./sum(n_b0[spcid,:],2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] = 1;
        n_fill = (prop_n .>= n_thresh)[:,1];
        #Build a list of those that pass each individual test
        a_pass = spcid[a_fill];
        n_pass = spcid[n_fill];
        #which species pass?
        survivors = intersect(a_pass,n_pass);
        spext1 = setdiff(spcid,survivors);
        
        #2) By not having the max strength for consuming at least one resource
        prext_comp = Array{Bool}(length(spcid));
        strength = sum(n_b0[spcid,[1;cid]],2) .- sum(a_b[spcid,[1;cid]],2);
        smatrix = Array{Float64}(copy(a_b[spcid,[1;cid]]));
        for i=1:length(spcid)
            smatrix[i,:] = a_b[spcid[i],[1;cid]] * strength[i];
        end
        smatrix[find(iszero,a_b[spcid,[1;cid]])] = NaN;
        for i=1:length(spcid)
            cmatrix = smatrix[:,find(!iszero,a_b[spcid[i],[1;cid]])];
            #Proportion of strength-max foods
            propmax = sum(cmatrix[i,:] .>= findmax(cmatrix,1)[1]')/length(cmatrix[i,:]);
            if propmax > 0
                prext_comp[i] = false;
            else
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
    end #end time steps
    R"plot($sprich,pch='.')"

    
    
end
