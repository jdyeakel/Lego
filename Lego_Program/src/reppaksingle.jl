function reppaksingle(S,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight)

  #Shared variables
  sprich = Array{Int64}(tmax,reps);
  rich = Array{Int64}(tmax,reps);
  conn = Array{Float64}(tmax,reps);
  ext_prim = Array{Int64}(tmax,reps);
  ext_sec = Array{Int64}(tmax,reps);
  pot_col = Array{Int64}(tmax,reps);
  cumid = Array{Int64}(tmax,reps);
  cumspid = Array{Int64}(tmax,reps);
  num_sp = Array{Int64}(reps);
  tfinal = Array{Int64}(reps);
  
  int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_species(S,probs,ppweight);
  
  num_play = length(diag(int_m));
  
  comgen = SharedArray{Int64}(reps,num_play,tmax);
  
  for r=1:reps
    
    num_sp[r] = length(find(x->x=='n',diag(int_m)));
    
    #Establish community template
    cid, c_m, crev_m, com_tp, com_tind, com_mp, com_mind = initiate_comm_func(int_m,tp_m,tind_m,mp_m,mind_m);
    status = "open";
    
    comid = (Array{Int64,1})[];
    push!(comid,cid);
    #rich = Array{Int64}(1);
    #rich[r,1] = length(comid);
    #sprich[r,1] = 1;
    
    t=0;
    while status == "open" && t <= tmax
      t = t+1;
      #Colonize with some probability
      rcol = rand();
      if rcol < rate_col && status == "open"
        status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,potcol = colonize_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
        #Update community evolution
        push!(comid,cid);
        #Update count of richness
        # push!(rich,length(cid));
      end
      
      sprich[t,r] = length(spcid);
      # length(unique(cid))-length(cid)
      conn[t,r] = (sum(com_tp))/(sprich[t,r]^2);
      rich[t,r] = length(cid);
      
      
      
      cumid[t,r] = sum(cid);
      cumspid[t,r] = sum(spcid);
      
      pot_col[t,r] = length(potcol);
      
      comgen[r,cid,t] = 1;
      
      
      
    end #end while loop
    
    tfinal[r] = t;
    #Record data important for understanding priority effects
    
  end #end repetition loop
  
  
  return(
  int_m,
  sprich,
  rich,
  conn,
  comgen,
  ext_prim,
  ext_sec,
  pot_col,
  num_sp,
  tfinal
  )
  
end
