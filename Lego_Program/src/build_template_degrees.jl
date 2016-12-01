function build_template_degrees(num_play, probs, ppweight)

  p_n=copy(probs[1]);
  p_a=copy(probs[2]);
  p_m=copy(probs[3]);
  p_i=copy(probs[4]); #Ignore with 1 - pr(sum(other))
  #Defining paiwise probabilities
  #These are probabilities of pairwise interactions within the whole universe of possible interactions (15)
  #THIS IS PROBABLY RIGHT
  pw_prob_init = [
    pr_na = p_n*(p_a/(p_a+p_n+p_i+p_m)) + p_a*(p_n/(p_a+p_i+p_n)),
    pr_nn = p_n*(p_n/(p_a+p_n+p_i+p_m)),
    pr_ni = p_n*(p_i/(p_a+p_n+p_i+p_m)) + p_i*(p_n/(p_a+p_n+p_i)),
    pr_nm = p_n*(p_m/(p_a+p_n+p_i+p_m)) + p_m*1, #(p_n/p_n),
    pr_ia = p_i*(p_a/(p_a+p_n+p_i)) + p_a*(p_i/(p_a+p_i+p_n)),
    pr_ii = p_i*(p_i/(p_a+p_n+p_i)),
    pr_aa = p_a*(p_a/(p_i+p_n+p_a))
  ];

  #Proposed method:
  #Build a random food web from the niche model given the number of species
  #Determine link-density properties from the web
  #Populate the interaction matrix with these statistical properties, but the structure will be scrambled. (I think)


  #Step 1: niche values
  nichev = rand(num_play);
  #Derive connectance from the pr('assimilate')*number of potential links
  #But first try to discount the number of players that are NOT species, which should be on average num_play*pr_nm
  S = num_play - (pr_nm*num_play);
  #Pr('A')*Directed Connectance
  C = ((pr_ia+pr_na+pr_aa)*(S*(S-1)))/(S^2);
  Ebeta = 2*C;
  beta = (1/Ebeta) - 1;
  BetaD = Beta(1,beta);
  rBeta = rand(BetaD,num_play);
  #NOTE: Perhaps I should accurately count how many prey fall within the range
  # rangev = length(find(x->x>nmin&&x<nmax,nichev));
  rangev = rBeta .* nichev;
  centerv = zeros(num_play);
  nmin = zeros(num_play);
  nmax = zeros(num_play);
  degrees = zeros(Int64,num_play);
  for i=1:num_play
    ud = Uniform(rangev[i]/2,nichev[i]);
    centerv[i] = rand(ud);
    nmin[i] = centerv[i] - (rangev[i]/2);
    nmax[i] = centerv[i] + (rangev[i]/2);
    degrees[i] = length(find(x->x>nmin[i]&&x<nmax[i],nichev));
    #Force everything to eat SOMETHING
    if degrees[i] == 0
      degrees[i] = 1;
    end
  end


  #Create an empty character array with dimensions equal to the number of players
  int_m = Array(Char,num_play,num_play);
  #Matrix of trophic-only interactions
  #The number of species in the trophic and mutualistic matrices will all start out the sames size as int_m, however, the objects will be trimmed later
  t_m = zeros(Int64,num_play,num_play);
  #Matrix of mutualistic interactions
  m_m = zeros(Int64,num_play,num_play);

  #Set array equal to zero
  int_m[1:(num_play*num_play)] = '0';

  prim_prod = zeros(num_play);
  #The first true species (row/col 2) is always a primary producer
  prim_prod[2] = 1;


  #Assigning trophic interactions to species in the net
  for i = 2:num_play

    #Assigning initial trophic interactions
    #We could add an additional vector of 1s to weight towards the basal resource
    #if we use a weighted vector, say 1/3 of species will be primary producers
    pvec = Array{Int64}(Int64(round(num_play*(ppweight))))*0+1;
    vec = collect(1:num_play);
    #Remove ith element
    deleteat!(vec,i);
    #weighted vec
    wvec = [pvec;vec];
    #What is the degree of player i?
    k = degrees[i]
    #Randomly choose the identities of prey without replacement
    if prim_prod[i] == 0
      resource = sample(wvec,k,replace=false)
    else
      #If you are a primary producer, at least 1 thing you conusume must be row/column 1
      resource_wop = sample(wvec,k-1,replace=false)
      resource = [resource_wop;1]
    end


    #Establish these prey in the interaction matrix
    int_m[i,collect(resource)] = 'a'
    #Initially, the resources are all species. Those that get turned into objects after 'm' interactions are deteremined will be trimmed from the trophic matrix (because only interactions with species are allowed)
    t_m[i,collect(resource)] = 1;
    #t_m[collect(resource),i] = 0;

    #Note: to vectorize a row from int_m (so that it is an Array{Char,1}), we would write v = int_m[i,:][:]

    #Assigning complimentary interactions
    #First disconsider aa interactions
    aa_int = find(x -> x == 'a', int_m[collect(resource),i])
    aa_int_loc = resource[aa_int];


    if length(aa_int) >= 1
      #Remove elements ee_int
      deleteat!(resource,aa_int)
    end

    #An 'a' has been drawn...
    #Now randomly determine ai or an (aa already emerges from assigning a's)
    #Draw from a binomial: which interaciton will be Asymmetric predation?
    #Draw 1 = Asymmetric predation; Draw 0 = facultative mutualism
    pr_i_given_a = (p_i/(p_i+p_n)); #again discounting a-a, which are already determined
    bindist = Binomial(1,pr_i_given_a)
    bin_draw = rand(bindist,length(resource))

    #Asymmetric predation
    res_i = resource[find(x->x==1,bin_draw)]
    int_m[collect(res_i),i] = 'i' #assigning n's according to random draw

    #Faculatative mutualism
    res_n = resource[find(x->x==0,bin_draw)]
    int_m[collect(res_n),i] = 'n' #assigning i's according to random draw

    #Record mutualistic interaction in m_m matrix
    m_m[collect(res_n),i] = 1;
    m_m[i,collect(res_n)] = 1;

  end

  #Duplicate this information for tall_m
  #This matrix includes indirect trophic interactions
  #Where if 1 (m) 2, and 3 (a) 2, then 3 (a) 1
  tall_m = copy(t_m);


  #Non-trophic interactions

  int_labels = ["na","nn","ni","nm","ia","ii","aa"]

  pw_prob_new = copy(pw_prob_init)
  #Eliminate n-a, i-a, a-a interactions, which are already determined
  deleteat!(pw_prob_new,[1,5,7])
  pw_prob_new = pw_prob_new/sum(pw_prob_new)

  new_int_labels = copy(int_labels)
  deleteat!(new_int_labels,[1,5,7])



  prob_line = cumsum(pw_prob_new);



  for i = 2:num_play
    for j = 2:num_play
      #Only choose interactions for empty elements (int_m = '0')
      #Only do this for the lower triangular part of the matrix
      if int_m[i,j] == '0' && i > j

        r_draw = rand()

        #N:N - symmetric mutualism
        if r_draw < prob_line[1]
          int_m[i,j] = 'n'
          int_m[j,i] = 'n'
          #Update mutualistic network
          m_m[i,j] = 1;
          m_m[j,i] = 1;
        end

        #N:I - commensalism
        if r_draw > prob_line[1] && r_draw < prob_line[2]
          index = 2 #find(x -> x == "ni",new_int_labels)
          mij = Array(Char,2)
          mij[1] = rand([new_int_labels[index][1],new_int_labels[index][2]])
          if mij[1] == new_int_labels[index][1]
            mij[2] = new_int_labels[index][2]
          else mij[2] = new_int_labels[index][1]
          end
          int_m[i,j] = mij[1]
          int_m[j,i] = mij[2]
        end


        #N:M - ecosystem engineering
        if r_draw > prob_line[2] && r_draw < prob_line[3]
          index = 3 #find(x -> x == "nm",new_int_labels)
          mij = Array(Char,2)
          mij[1] = rand([new_int_labels[index][1],new_int_labels[index][2]])
          if mij[1] == new_int_labels[index][1]
            mij[2] = new_int_labels[index][2]
          else mij[2] = new_int_labels[index][1]
          end
          int_m[i,j] = mij[1]
          int_m[j,i] = mij[2]
          #Ensure nothing 'makes' the forced primary producer
          #Not sure if I need this 10/24/16
          if j == 2
            int_m[i,j] = 'n'
            int_m[j,i] = 'm'
          end
        end


        #I:I - neutral interaction
        if r_draw > prob_line[3] && r_draw < prob_line[4]
          int_m[i,j] = 'i'
          int_m[j,i] = 'i'
        end

      end #if
    end #j
  end #i

  #Filling in the diagonal
  #Determine which are active players and which aren't
  for i = 1:num_play
    if length(find(x->x=='a',int_m[i,:])) > 0
      int_m[i,i] = 'n'
    else int_m[i,i] = 'i'
    end
  end


  #Implement posthoc Rules
  #1: Row/Col 1 is the sun
  int_m[1,:] = 'i'
  #All column elements that ARENT 'a' are 'ignore'
  force_ignore=find(x->x!='a',int_m[:,1]);
  int_m[force_ignore,1] = 'i';

  #due to the above forcing, there are no mutualistic interactions with the sun
  #The only 'a' interactions with the sun are 'a-i', not 'a-n'
  m_m[1,:] = 0; #this probably isn't needed.
  m_m[:,1] = 0;

  #A matrix for Direct + Indirect mutualistic interactions
  mall_m = copy(m_m);


  #Basal primary producer (Row/Col 2) doesn't need anything
  colonizer_n=find(x->x=='n',int_m[2,:]);
  int_m[2,colonizer_n[2:length(colonizer_n)]] = 'i';

  #Connectance before excluding 'made' things
  # Sorig = length(find(x->x=='n',diag(int_m)));
  # Lorig = sum(t_m);
  # Corig = Lorig/(Sorig^2);

  #2: if player A contains an "m" with player B, player B is "i" with everything
  # except if it is "n" with player A

  #NOTE: multiple A,B,C can make the same D
  #NOTE: Should 'made' things only be needed and not assimilated, which implies biomass flow? i.e. trophic int?
  
  #BUG: There are 'objects' that need other objects...
  # REASON: if a high R/C species (10) makes an object lower on the list (5), and a species intermediate (8) makes (10), then, object 5 needs species 10, but then species 10 get turned into an object... so the object needs an object!
  
  madev = Array{Int64}(0);
  for i = 2:num_play
    int_v = copy(int_m[i,:]);
    made = find(x->x=='m',int_v);
    l_made = length(made);
    if l_made > 0
      for k = 1:l_made
        #Compile list of made things
        push!(madev,made[k])
        
        #If the made thing 'makes' anything, the reverse 'need' interactions must be turned to 'ignore'
        unmade = find(x->x=='m',int_m[made[k],:]);
        int_m[unmade,made[k]] = 'i'
        
        #The made thing ignores everything... including diag
        int_m[made[k],:] = 'i'

        #what thing(s) make it? ~ Restablishes the m-n interaction
        makers = find(x->x=='m',int_m[:,made[k]]);
        int_m[made[k],makers] = 'n'; #Except the thing(s) that makes it

        #For the trophic matrices, which include only species interactions, the made object interacts with nothing except the species that makes it
        #It is not a living thing, so nothing 'eats' it
        t_m[made[k],:] = 0;
        t_m[:,made[k]] = 0;
        tall_m[made[k],:] = 0;
        tall_m[:,made[k]] = 0;
        m_m[made[k],:] = 0;
        m_m[:,made[k]] = 0;

      end
    end
  end

  #Documenting indirect trophic and mutualistic interactions
  #We want to do this outside of the above maker loop because we want all species/objects correctly accounted for before we try to count the indirect interactions

  obrc = find(x->x=='i',diag(int_m));
  for i=2:length(obrc) #Start at 2 because nothing makes the basal resource
    
    #Define the made object
    made = obrc[i]
    #Define all species that make the object
    makers = find(x->x=='m',int_m[:,made]);
    
    #Indirect trophic interactions
    
    #NOTE: If a species eats a made object, it is indirectly 'eating' the species that made the object. We count this as an indirect interaction and record it in tall_m, and this gets trimmed to tind_m.
    #The consumers of the thing that species i makes indirectly consume species i
    #What species EAT the made things?
    ind_cons = find(x->x=='a',int_m[:,made]);
    tall_m[ind_cons,makers] = 1;

    #Indirect mutualistic interactions

    #We have to check
    #1) If species 2 eat/needs an object that species 1 makes
    #2) If species 1 needs/eats species 2
    #I would suspect that this would be rare.
    ind_need = find(x->x=='n'||x=='a',int_m[:,made]);
    for j=1:length(ind_need)
      #Cycle through each thing that makes object to assess different indirect interactions
      for k=1:length(makers)
        
        #If Species 2 (ind_need[j]) NEEDS Object (made)
        #AND Speces 1 (makers[k]) EATS Species 2 (ind_need[j])
        #Then we establish a mutualistic interaction between Species 1 (makers [k]) and Species 2 (ind_need[j])
        if int_m[ind_need[j],made] == 'n' && int_m[makers[k],ind_need[j]] == 'a' 
          mall_m[ind_need[j],makers[k]] = 1;
          mall_m[makers[k],ind_need[j]] = 1;
        end
        
        #If Species 2 (ind_need[j]) EATS Object (made)
        #AND Speces 1 (makers[k]) NEEDS Species 2 (ind_need[j])
        #Then we establish a mutualistic interaction between Species 1 (makers [k]) and Species 2 (ind_need[j])
        if int_m[ind_need[j],made] == 'a' && int_m[makers[k],ind_need[j]] == 'n' 
          mall_m[ind_need[j],makers[k]] = 1;
          mall_m[makers[k],ind_need[j]] = 1;
        end
        
        #If Species 2 (ind_need[j]) NEEDS Object (made)
        #AND Speces 1 (makers[k]) NEEDS Species 2 (ind_need[j])
        #Then we establish a mutualistic interaction between Species 1 (makers [k]) and Species 2 (ind_need[j])
        if int_m[ind_need[j],made] == 'n' && int_m[makers[k],ind_need[j]] == 'n' 
          mall_m[ind_need[j],makers[k]] = 1;
          mall_m[makers[k],ind_need[j]] = 1;
        end
        
      end
    end
  end

  # #MANY I-I rows/columns? nope
  # itestrc = Array{Bool}(num_play);
  # itestr = Array{Bool}(num_play);
  # itestc = Array{Bool}(num_play);
  # for i=1:num_play
  #   itestr[i] = all(int_m[i,:].=='i');
  #   itestc[i] = all(int_m[:,i].=='i');
  #   itestrc[i] = (itestr[i]== true && itestc[i]==true);
  # end

  #Which rows/columns are species?
  sprc = find(x->x=='n',diag(int_m));
  num_sp=length(sprc);
  #Build empty matrices to store species-only direct/indirect trophic/mutualistic interactions
  #The +1 is due to adding Row/Col 1 for the basal resource
  sp_m = Array{Char}(num_sp+1,num_sp+1);
  tp_m = zeros(Int64,num_sp+1,num_sp+1);
  tind_m = zeros(Int64,num_sp+1,num_sp+1);
  mp_m = zeros(Int64,num_sp+1,num_sp+1);
  mind_m = zeros(Int64,num_sp+1,num_sp+1);
  #ensure the 1st row/column in species interaction matrices is the sun (basal resource)
  #This means the 1st row is 'i' for sp_m, or 0 for tp_m, tind_m, mp_m
  #This means the 1st column is

  #NOTE: I'm not sure why, but this vec() is important. It's an Array splicing issue I think.
  sp_m[1,:] = ['i';int_m[1,sprc]];
  sp_m[:,1] = ['i';int_m[sprc,1]];
  tp_m[1,:] = [0;vec(t_m[1,sprc])];
  tp_m[:,1] = [0;vec(t_m[sprc,1])];
  tind_m[1,:] = [0;vec(tall_m[1,sprc])];
  tind_m[:,1] = [0;vec(tall_m[sprc,1])];
  mp_m[1,:] = [0;vec(m_m[1,sprc])];
  mp_m[:,1] = [0;vec(m_m[sprc,1])];
  mind_m[1,:] = [0;vec(mall_m[1,sprc])];
  mind_m[:,1] = [0;vec(mall_m[sprc,1])];


  #The rest of the matrices (2:num_sp+1 - the +1 is due to the fact that we have placed the basal resource in row/column 1) will be filled in with the species only part of the full matrices
  sp_m[collect(2:num_sp+1),collect(2:num_sp+1)] = copy(int_m[sprc,sprc]);
  tp_m[collect(2:num_sp+1),collect(2:num_sp+1)] = copy(t_m[sprc,sprc]);
  tind_m[collect(2:num_sp+1),collect(2:num_sp+1)] = copy(tall_m[sprc,sprc]);
  mp_m[collect(2:num_sp+1),collect(2:num_sp+1)] = copy(m_m[sprc,sprc]);
  mind_m[collect(2:num_sp+1),collect(2:num_sp+1)] = copy(mall_m[sprc,sprc]);

  # #How many 'made things?'
  # Snew = length(find(x->x=='n',diag(int_m)));
  # Lnew = sum(t_m)/2;
  # Cnew = Lnew/(Snew^2);


  #NOTE: This takes a long time
  #NOTE: if par is true, it will calculate similarity in parallel; if 2000 species/objects are in the matrix, it shaves off 100 seconds
  #NOTE: I could make this much faster if I only calculate similarity between species (not objects). I'm currently doing it for both. This could cut it by half.
  
  #Optional code to build similarity matrix
  #simvalue = zeros(num_play,num_play);
  
  #We need simvalue to be declared before the if loop.
  #It is redefined as an Array or SharedArray in the loop
  
  simvalue = Array{Float64}(num_play,num_play);
  # simvalue = 0;
  # 
  # if par == true
  #   simvalue = SharedArray{Float64}(num_play,num_play);
  #   if sim == true
  #     @sync @parallel for i = 1:num_sp
  #       for j = 1:num_sp
  #         if j > i
  #           seq1 = copy(int_m[sprc[i],:]);
  #           seq2 = copy(int_m[sprc[j],:]);
  #           similarity = sim_func(seq1,seq2,num_play);
  #           #The 0.00001 is to ensure a nonzero value for any species
  #           simvalue[sprc[i],sprc[j]] = copy(similarity)+0.00001;
  #           simvalue[sprc[j],sprc[i]] = copy(similarity)+0.00001;
  #         end
  #         if j == i
  #           simvalue[sprc[i],sprc[j]] = 1.0;
  #         end
  #       end
  #     end
  #   end
  # else
  #   simvalue = zeros(num_play,num_play);
  #   if sim == true
  #     for i = 1:num_sp
  #       for j = 1:num_sp
  #         if j > i
  #           seq1 = copy(int_m[sprc[i],:]);
  #           seq2 = copy(int_m[sprc[j],:]);
  #           similarity = sim_func(seq1,seq2,num_play);
  #           #The 0.00001 is to ensure a nonzero value for any species
  #           simvalue[sprc[i],sprc[j]] = similarity+0.00001;
  #           simvalue[sprc[j],sprc[i]] = similarity+0.00001;
  #         end
  #         if j == i
  #           simvalue[sprc[i],sprc[j]] = 1.0;
  #         end
  #       end
  #     end
  #   end
  # end
  # 


return(int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue)

end
