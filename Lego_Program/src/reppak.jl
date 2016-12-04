function repsim(num_play,reps,tmax,a_thresh,n_thresh,trophicload,rate_col,probs,ppweight,sim,par)

  #Shared variables
  sprich = SharedArray{Int64}(tmax,reps);
  rich = SharedArray{Int64}(tmax,reps);
  conn = SharedArray{Float64}(tmax,reps);
  comgen = SharedArray{Int64}(reps,num_play,tmax);
  ext_prim = SharedArray{Int64}(tmax,reps);
  ext_sec = SharedArray{Int64}(tmax,reps);
  pot_col = SharedArray{Int64}(tmax,reps);
  
  @time int_m, sp_m, t_m, tp_m, tind_m, mp_m, mind_m, simvalue = build_template_degrees(num_play,probs,ppweight,sim,par);
  
  @sync @parallel for r=1:reps
    #Establish community template
    cid, c_m, crev_m, com_tp, com_tind, com_mp, com_mind = initiate_comm_func(int_m,tp_m,tind_m,mp_m,mind_m);
    status = "open";
    
    comid = (Array{Int64,1})[];
    push!(comid,cid);
    rich = Array{Int64}(1);
    rich[1] = length(comid);
    
    while status = "open"
      
      #Colonize with some probability
      rcol = rand();
      if rcol < rate_col && status == "open"
        status,cid,spcid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind,potcol = colonize_func(int_m,tp_m,tind_m,mp_m,mind_m,a_thresh,n_thresh,cid,c_m,crev_m,com_tp,com_tind,com_mp,com_mind);
        #Update community evolution
        push!(comid,cid);
        #Update count of richness
        push!(rich,length(cid));
      end
      
      sprich[t,r] = length(spcid);
      # length(unique(cid))-length(cid)
      conn[t,r] = (sum(com_tp))/(sprich[t,r]^2);
      rich[t,r] = length(cid);
      
      pot_col[t,r] = length(potcol);
      
      comgen[r,cid,t] = 1;
      
      
    end #end while loop
    
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
  pot_col
  )
  
end
