function colext(
    int_m,tp_m,tind_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,
    cid,a_thresh,n_thresh,extmid,steep,epsilon,sigma,prext,colcheck,extcheck,exttype)
    
    
    ######################
    #CURRENT COMMUNITY
    ######################
    
    #Define the status of the community as open to colonizers
    status = 1;

    #ID's of agents in the community
    cid_old = copy(cid);

    #IDs of species only
    spcid = intersect(sp_v,cid_old);
    
    colrand = rand();
    if colcheck > colrand
    
        ######################
        #COLONIZATION MODULE
        ######################

        #which species have an 'a'-link to anyone in the community? (this should include primary producers and those that are trophically linded to objects: [1 cid])
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
        pot_col = intersect(a_pass,n_pass);

        #Count the number that pass
        lpot_col = length(pot_col);

        #Choose a colonizer at random
        if lpot_col > 0
            c_sp = rand(pot_col);

            #select made objects as well
            c_ob_full = find(isodd,m_b[c_sp,:]);

            #Only bring in objects that aren't already in the community
            c_ob = setdiff(c_ob_full,cid_old);

            colonize = [c_sp;c_ob];

            #Update community
            cid = [cid_old;colonize];
        else
            status = 0; #status=0 means that the community is closed
        end
    else #If colonization is closed
        lpot_col = 0;
        cid = copy(cid_old);
    end #End colonization module



    ######################
    #EXTINCTION MODULE
    ######################
    num_ext1 = 0;
    num_ext2 = 0;

    #run the extinction module or colonization by itself?
    extrand = rand();
    if extcheck > extrand

        cid_old = copy(cid);

        #which are species?
        spcid = intersect(sp_v,cid_old);
        #which are objects?
        ocid = setdiff(cid_old,spcid);

        #1) determine predation load for each species
        # num_preds = sum(a_b[spcid,spcid],1)[1,:];

        #2) Caluculate extinction probability
        # avgk = 8; # avgk = gL/gS;
        # prext_pred = (1./(1+exp.(-0.5.*(num_preds.-(trophicload*avgk)))));
        
        #Based on the Holt 1977 model where N* propto 1/n where n is the number of preds
        #baseline extinction probability
        
        binext = Array{Int64}(length(spcid));
        num_ext1 = Int64;
        
        #Number of predators consuming each species
        num_preds = sum(a_b[spcid,spcid],1)[1,:];

        
        if exttype == "PL"
            #Extinction by predation load
            #=========================#
            # epsilon = 0.01;
            # sigma = 1/3;
            prext_pred = 0.5*erfc.((1-num_preds*epsilon)/(sigma*sqrt(2)));
            #3) Draw extinctions
            binext = rand.(Binomial.(1,prext_pred));
            num_ext1 = sum(binext);
            #=========================#

        end
        if exttype == "RO"
            
            #Extinction by similarity
            #=========================#
            
            #Species consuming each resource (species and objects)
            # preds = vec(sum(a_b[cid,cid],1));
            #Species consuming or needing each resource (species and objects)
            users = vec(sum(a_b[cid,cid],1)) .+ vec(sum(n_b0[cid,cid],1));
            #Number of resources for each species (will be zero if species are only eating basal resource)
            # res = vec(sum(a_b[spcid,cid],2));
            used = vec(sum(a_b[spcid,cid],2)) .+ vec(sum(n_b0[spcid,cid],2));
            
            #Proporitonal overlap of resources between species
            # res_overlap = (((a_b[spcid,cid]*preds).-res)/(length(spcid)-1))./res;
            res_overlap = (((((a_b[spcid,cid]+n_b0[spcid,cid]))*users).-used)/(length(spcid)-1))./used;
            #nan if 0 resources ~ pure primary producer
            res_overlap[isnan.(res_overlap)] = 0;
            
            #Similarity where pr(extinction) = 0.5
            extmidpoint = extmid;
            abeta = steep; #Higher values = steeper sigmoid
            
            bbeta = (-1 + 3*abeta + 2*extmidpoint - 3*abeta*extmidpoint)/(3*extmidpoint);
            pr_background = 0.000;
            # R"plot($(collect(0.0:0.001:1.0)),$(pr_background + (1-pr_background)*cdf.(Beta(abeta,bbeta),collect(0.0:0.001:1.0))),pch='.',log='x'); points($extmidpoint,0.5,pch=16)"
            prext_comp = pr_background + (1-pr_background)*cdf.(Beta(abeta,bbeta),res_overlap);
            binext = rand.(Binomial.(1,prext_comp));
            num_ext1 = sum(binext);
            #=========================#
            
            # if exttype == "ROPL"
            #     #TODO: Merge the RO extinction metric to the predation load extinction metric
            # end
            
        end
        if exttype == "EN"
            prext_comp = Array{Float64}(length(spcid));
            #Calculate the strength of the interaction (num needs - num eats)
            strength = sum(n_b0[spcid,[1;cid]],2) .- sum(a_b[spcid,[1;cid]],2);
            smatrix = Array{Float64}(copy(a_b[spcid,[1;cid]]));
            for i=1:length(spcid)
                smatrix[i,:] = a_b[spcid[i],[1;cid]] * strength[i];
            end
            smatrix[find(iszero,a_b[spcid,[1;cid]])] = NaN;
            for i=1:length(spcid)
                
                cmatrix = smatrix[:,find(!iszero,a_b[spcid[i],[1;cid]])]
                
                #Proportion of strength-max foods
                propmax = sum(cmatrix[i,:] .>= findmax(cmatrix,1)[1]')/length(cmatrix[i,:]);
                
                if propmax > 0
                    prext_comp[i] = 0;
                else
                    prext_comp[i] = prext;
                end
                #probability of extinction is 1 - proportion of strength-max foods
                # prext_comp[i] = 1 - propmax;
                
                # if propmax > prext
                #     prext_comp[i] = prext;
                # else
                #     prext_comp[i] = 1 - prext;
                # end
            end
            binext = rand.(Binomial.(1,prext_comp));
            num_ext1 = sum(binext);
        end

        #Only go through this round if there are any extinctions
        if num_ext1 > 0

            #####################
            # Primary extinctions
            #####################

            #Which species are going extinct?
            e_sp = spcid[find(isodd,binext)];
            le_sp = length(e_sp);

            #Delete species and their UNIQUELY MADE objects
            #what objects are uniquely made by esp?
            # e_ob = ocid[(sum(m_b[e_sp,ocid],1).==1)[1,:]];
            # e_ob = ocid[find(x->x==1,sum(m_b[e_sp,ocid],1))];
            #The minimum will be 1... subtract 1 from all and use find(iszero) for speed
            # e_ob = ocid[find(iszero,sum(m_b[e_sp,ocid],1)-1)];
            e_ob = ocid[find(iszero,sum(m_b[spcid,ocid],1)-sum(m_b[e_sp,ocid],1))]
            le_ob = length(e_ob);

            prim_extinct = [e_sp;e_ob];

            #Update community
            cid = setdiff(cid_old,prim_extinct);
            
            #which are species?
            spcid = intersect(sp_v,cid);
            #which are objects?
            ocid = setdiff(cid,spcid);
            spcid_ind = indexin(spcid,[1;sp_v]);


            #######################
            # Secondary extinctions
            #######################

            #Only do this if there are still members in the community
            if length(cid) > 0
                
                #Determine if the SUN has been disconnected
                
                adjmatrix = tind_m[[1;spcid_ind],[1;spcid_ind]];
                
                
                # #====================#
                # #Connected component method
                # #Create a symmetric matrix of the adjacency
                # symmatrix = (adjmatrix + adjmatrix').>0;
                # gs = SimpleGraph(symmatrix);
                # #Find the connected components and eliminate anything that is not a component with vertex 1;
                # cc = connected_components(gs);
                # #Collect the disconnected species
                # todisconnect = cc[find(x->in(1,x)==false,cc)];
                
                #====================#
                #Path to sun method
                g = DiGraph(adjmatrix');
                #Find directed pathlengths to sun
                # paths = deleteat!(length.(enumerate_paths(dijkstra_shortest_paths(g,1))),1)-1;
                #FASTEST algorithm I've found
                paths = gdistances(g,1);
                #Basal resource itself will have a length of '0'
                #Disconnected nodes have an insanely large path...
                #If any nodes to NOT have a directed path length, then they must be disconnected.
                todisconnect = find(x->x>size(adjmatrix)[1],paths);
                
                
                conn_test = true;
                if length(todisconnect) > 0
                    disconnected = collect(Iterators.flatten(todisconnect));
                    
                    #convert these identifiers (elements of tind_m[[1;spcid_ind],[1;spcid_ind]]) back to species IDs
                    disconnected = [1;spcid_ind][disconnected];
                    
                    conn_test = false;
                else
                    disconnected = Array{Int64}(0);
                end
                
                # basalconsumer = sum(a_b[[1;cid],[1;cid]][:,1]);
                
                if conn_test == false
                    
                    survivors = setdiff(spcid,disconnected);
                    e_sp2 = disconnected;

                    #Number of secondary extinctions
                    le_sp2 = length(e_sp2);
                    num_ext2 += le_sp2;
                    
                    e_ob2 = ocid[find(iszero,sum(m_b[spcid,ocid],1)-sum(m_b[e_sp2,ocid],1))];

                    le_ob2 = length(e_ob2);

                    sec_extinct = [e_sp2;e_ob2];

                    #Update community
                    cid = setdiff(cid,sec_extinct);
                    
                else

                    #The assimilate/need test has to pass (an_test == true) to avoid any more secondary extinctions
                    an_test = false;

                    while an_test == false
                        #Which are species?
                        spcid = intersect(sp_v,cid);
                        #Which are objects?
                        ocid = setdiff(cid,spcid);

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

                        #If there are secondary extinctions...
                        if length(survivors) != length(spcid)

                            #Species that secondarily don't pass the a-n test
                            e_sp2 = setdiff(spcid,survivors);

                            #Number of secondary extinctions
                            le_sp2 = length(e_sp2);

                            num_ext2 += le_sp2;

                            #Delete species and their UNIQUELY MADE objects
                            #what objects are uniquely made by esp?
                            # e_ob2 = ocid[(sum(m_b[e_sp2,ocid],1).==1)[1,:]];
                            # e_ob2 = ocid[find(x->x==1,sum(m_b[e_sp2,ocid],1))];

                            # e_ob2 = ocid[find(iszero,sum(m_b[e_sp2,ocid],1)-1)];
                            # sp_ob = ocid[find(iszero,sum(m_b[spcid,ocid],1)-1)];


                            e_ob2 = ocid[find(iszero,sum(m_b[spcid,ocid],1)-sum(m_b[e_sp2,ocid],1))];

                            le_ob2 = length(e_ob2);

                            sec_extinct = [e_sp2;e_ob2];

                            #Update community
                            cid = setdiff(cid,sec_extinct);

                        else

                            #If there are no survivors, test is passed; exit loop
                            an_test = true;

                        end

                    end #end secondary extinction while loop
                
                end #end sun-disconnected loop 

            end #end if length(cid) > 0

        end #end extinctions module

    end #end (if extinctions==true) module

    return(
    cid,
    lpot_col,
    status,
    num_ext1,
    num_ext2
    )

end #end colext function
