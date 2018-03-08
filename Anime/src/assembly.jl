function assembly(a_thresh,n_thresh,trophicload,extinctions,tmax)
    
    cid = Array{Int64}(0);
    rich = Array{Int64}(tmax);
    sprich = Array{Int64}(tmax);
    turnover = Array{Float64}(tmax);
    prim_ext = Array{Int64}(tmax);
    sec_ext = Array{Int64}(tmax);
    
    for t = 1:tmax
      #Print every 1000 timesteps
      if mod(t,1000)==0
        println(string("t=",t))
      end
      
      cid_old = copy(cid);
      
      cid,
      lpot_col,
      status,
      num_ext1,
      num_ext2 = colext(int_m,cid,a_thresh,n_thresh,extinctions,trophicload);
      
      rich[t], sprich[t], turnover[t] = dynstructure(cid,cid_old);
      
      prim_ext[t] = num_ext1;
      sec_ext[t] = num_ext2;

    end
    
    return(cid,rich,sprich,turnover,prim_ext,sec_ext)
    
end