#Does not generate comgen output, and does NOT calculate potential colonizers (faster)

function repsimintsingle(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight,runpsw)

  #Shared variables
  sprich = SharedArray{Int64}(tmax,reps);
  rich = SharedArray{Int64}(tmax,reps);
  conn = SharedArray{Float64}(tmax,reps);
  #comgen = SharedArray{Int64}(reps,num_play,tmax);
  ext_prim = SharedArray{Int64}(tmax,reps);
  ext_sec = SharedArray{Int64}(tmax,reps);
  tw = SharedArray{Float64}(tmax,reps);
  twind = SharedArray{Float64}(tmax,reps);
  #int_mv = SharedArray{Char}(num_play*reps,num_play);
  psw = SharedArray{Float64}(tmax,reps);
  pswind = SharedArray{Float64}(tmax,reps);

  @sync @parallel for r=1:reps

    int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_species(S,probs,ppweight);

    num_play = length(diag(int_m));


    #Establish community template
    cid, c_m, crev_m, com_tp, com_tind, com_mp, com_mind = initiate_comm_func(int_m,tp_m,tind_m,mp_m,mind_m);
    
    fwt = Array(Array{Int64},tmax);
    fwtind = Array(Array{Int64},tmax);
    fwm = Array(Array{Int64},tmax);
    fwmind = Array(Array{Int64},tmax);

    timetic=0;
    for t=1:tmax
      timetic = timetic+1;
      status = "open";
      #Colonize with some probability
      rcol = rand();
      if rcol < rate_col && status == "open"
        status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind = colonizesingle_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
      end
      #Always run extinction module because probabilities are assessed within
      status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,extinctions = extinct_func2(int_m,tp_m,a_thresh,n_thresh,trophicload,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
      # status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,extinctions = extinct_func(int_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,simvalue);
      #Save primary and secondary extinction information
      ext_prim[t,r] = extinctions[1];
      ext_sec[t,r] = extinctions[2];
      sprich[t,r] = length(spcid);
      # length(unique(cid))-length(cid)
      conn[t,r] = (sum(com_tp))/(sprich[t,r]^2);
      rich[t,r] = length(cid);
      tw[t,r] = mean(trophicwidth(com_tp));
      twind[t,r] = mean(trophicwidth(com_tind));
      
      if runpsw == true
          #Record interaction matrices
          fwt[t] = copy(com_tp);
          fwtind[t] = copy(com_tind);
          fwm[t] = copy(com_mp);
          fwmind[t] = copy(com_mind);
      end
      
      
      #comgen[r,cid,t] = 1;
    end #end time loop
    #int_mv[1+(r-1)*num_play:num_play*r,1:num_play] = int_m;
    
    if runpsw == true
        
        ##########################
        #Proportion of stable webs
        ##########################
        psw_r = zeros(Float64,tmax);
        pswind_r = zeros(Float64,tmax);
        for t=1:tmax
            sigma = 0.4
            reps = 50;
            psw_r[t] = PSWebs(fwt[t],fwm[t],sigma,reps);
            pswind_r[t] = PSWebs(fwtind[t],fwmind[t],sigma,reps);
            # if mod(t,100)==0
            #   println(string("t=",t))
            # end
        end
        psw[:,r] = psw_r;
        pswind[:,r] = pswind_r;
        
    end
    
    
  end #end repetition loop



  return(
  #int_mv,
  sprich,
  rich,
  conn,
  #comgen,
  ext_prim,
  ext_sec,
  tw,
  twind,
  psw,
  pswind
  )

end
