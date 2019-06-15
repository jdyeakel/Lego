function assembly_delayedobjects(int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,lambda,
    athresh,nthresh,maxits,cn,ce,cp)

    S = length(sp_v) + 1;
    N = size(int_m)[1];
    MaxN = convert(Int64,floor(S + S*lambda));
    cid = Array{Int64}(undef,0);
    sprich = Array{Int64}(undef,0);
    rich = Array{Int64}(undef,0);
    clock = Array{Float64}(undef,0);
    CID = falses(N,maxits);

    #NOTE strength matrix can't be built a-priori!
    # #Build the strength matrix apriori
    # strength = vec(pi*sum(n_b0,dims=2)) .- vec(sqrt(2)*sum(a_b,dims=2)) .- vec(sum(a_b,dims=1));
    # smatrix = Array{Float64}(copy(a_b));
    # for i=1:size(a_b)[1]
    #     smatrix[i,:] .= a_b[i,:] * strength[i];
    # end
    # # smatrix[findall(iszero,a_b)] = NaN;
    #
    minstrength = -ce*Float64(S) - cp*Float64(S);


    t=0;
    it = 0;
    while it < maxits

        #does sorting CID make a difference?
        sort!(cid);

        cid_old = copy(cid);
        #Which are species?
        spcid = intersect(sp_v,cid);
        #Which are objects?
        ocid = setdiff(cid,spcid);


        #COUNT POTENTIAL COLONIZERS
        trophiclinked = setdiff(int_id[(sum(a_b[:,[1;cid]],dims=2) .> 0)[:,1]],cid);
        #For each trophiclinked, count number of assimilate and need interactions in system
        #Determine in the proportion that already exists is >= the threshold
        a_fill = ((sum(a_b[trophiclinked,[1;cid]],dims=2)./sum(a_b[trophiclinked,:],dims=2)) .>= athresh)[:,1];
        prop_n = sum(n_b0[trophiclinked,[1;cid]],dims=2)./sum(n_b0[trophiclinked,:],dims=2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] .= 1;
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
        a_fill = ((sum(a_b[spcid,[1;cid]],dims=2)./sum(a_b[spcid,:],dims=2)) .> athresh)[:,1];
        prop_n = sum(n_b0[spcid,[1;cid]],dims=2)./sum(n_b0[spcid,:],dims=2);
        #If there are no 'need' interactions, proportion filled is always 1, and will always pass
        prop_n[isnan.(prop_n)] .= 1;
        n_fill = (prop_n .>= nthresh)[:,1];
        #Build a list of those that pass each individual test
        a_pass = spcid[a_fill];
        n_pass = spcid[n_fill];
        #which species pass?
        survivors = intersect(a_pass,n_pass);
        spext1 = setdiff(spcid,survivors);

        if length(spcid) > 0

            # Build the strength matrix at each community state
            # Strength values change over time so need to ve updated
            # Only record strength values of species (hence the [1:length(spcid)])
            #NOTE: Needs won't change; Eats is based on POTENTIAL niche; Vuln changes per timestep
            strength = vec(cn*sum(n_b0[spcid,cid],dims=2)) .- vec(ce*sum(a_b[spcid,:],dims=2)) .- (vec(cp*sum(a_b[spcid,cid],dims=1))[1:length(spcid)]);

            #This matrix applies the strength values across foraging interactions for each species in an spcid x spcid array
            #It does NOT include foraging interactions on the sun!
            cmatrix = zeros(Float64,length(spcid),length(cid));
            for i=1:length(spcid)
                cmatrix[i,:] .= a_b[spcid[i],cid] * strength[i];
            end

            # cmatrix = a_b[spcid,cid] .* reshape(repeat(strength,outer=length(cid)),length(spcid),length(cid));

            # cmatrix = (a_b[spcid,cid]' * (Matrix{Float64}(I,length(strength),length(strength)) .* strength))'

            #'zero' entrees need to be lower than any possible strength
            #So they will be effectively ignored
            #ce*S-cp*S is the theoretical min strength
            cmatrix[cmatrix.==0] .= minstrength;

            #Maximum strength values for each resource utilized
            cmax = findmax(cmatrix,dims=1)[1];

            prext_comp = trues(length(spcid));
            for i=1:length(spcid)

                #Don't count sun
                #catalogue prey for all species/objects in the system (excluding sun)
                ieats = Array{Bool}(a_b[spcid[i],cid]);
                #If you have >= the max strength for any of those prey, you stay
                #This means that a pure primary producer is not evaluated
                #any(ieats) = FALSE for pure primary producers;;;this means that primary producers do not go extinct by this method

                #For NON PRIMARY PRODUCERS
                #(any(strength[i] .>= cmax[ieats])==false) is TRUE if there is not any strength that is not >= the max strength on that food -> true prext
                #(any(strength[i] .>= cmax[ieats])==false) is FALSE if there is one or more strengths that is >= the max strength on that food -> false prext
                prext_comp[i] = any(ieats)*(any(strength[i] .>= cmax[ieats])==false);

            end
            #Grab those species destined for extinction
            spext2 = spcid[prext_comp];
        else
            spext2 = Array{Int64}(undef,0);
        end


        #
        # spext2 = Array{Int64}(0);
        spext = unique([spext1;spext2]);
        lspext = length(spext);


        #COUNT POTENTIAL EXTINCT OBJECTS
        obext = ocid[findall(iszero,vec(sum(m_b[spcid,ocid],dims=1)))];
        lobext = length(obext);
        levents = sum([lcol;lspext;lobext]);

        dt = 1/levents;

        t += dt;
        it += 1;
        push!(clock,t);
        
        #NOTE: Reduce probability that objects are chosen to be eliminated
        # skipobjectweight = 0.1;
        # lobext = Int64(floor(lobext*skipobjectweight));
        # levents = sum([lcol;lspext;lobext]);

        #Choose a random event
        re = rand();

        if re < (lcol/levents)

            #COLONIZATION FUNCTION
            c_sp = rand(col);
            #select made objects as well
            c_ob_full = findall(isodd,m_b[c_sp,:]);
            #Only bring in objects that aren't already in the community
            c_ob = setdiff(c_ob_full,cid_old);
            colonize = [c_sp;c_ob];
            #Update community
            cid = [cid_old;colonize];
            #Update species
            spcid = [spcid;c_sp];
        end

        if re > (lcol/levents) && re < ((lcol + lspext)/levents)

            #SPECIES EXTINCTION FUNCTION
            #select species to go extinct
            sp_bye = rand(spext,1);
            cid = setdiff(cid_old,sp_bye);
            spcid = setdiff(spcid,sp_bye);
        end

        if re > ((lcol + lspext)/levents)
            
            #What is the probability that it is eliminated?
            #This is to institute a 'delay' in object extinction
            #without biasing the Gillespie algorithm
            prob_elim = 0.5;
            if rand() < prob_elim
                #OBJECT EXTINCTION FUNCTION
                ob_bye = rand(obext,1);
                cid = setdiff(cid_old,ob_bye);
            end

        end

        #NOTE - updating CID....
        CID[cid,it] .= true;
        push!(sprich,length(spcid));
        push!(rich,length(cid));
    end #end time steps

    return(
    sprich,
    rich,
    clock,
    CID
    )

end


