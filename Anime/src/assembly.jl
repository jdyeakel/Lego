function assembly(
    int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tind_m,
    a_thresh,n_thresh,extinctions,tmax)
    
    cid = Array{Int64}(0);
    rich = Array{Int64}(tmax);
    sprich = Array{Int64}(tmax);
    turnover = Array{Float64}(tmax);
    prim_ext = Array{Int64}(tmax);
    sec_ext = Array{Int64}(tmax);
    res_overlap = Array{Float64}(tmax);
    conn = Array{Float64}(tmax);
    status = Array{Int64}(tmax);
    lpot_col = Array{Int64}(tmax);
    
    for t = 1:tmax
      # if mod(t,1000)==0
      #   println(string("t=",t))
      # end
      
      cid_old = copy(cid);
      
      cid,
      lpot_col[t],
      status[t],
      prim_ext[t],
      sec_ext[t] = colext(
      int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,
      cid,a_thresh,n_thresh,extinctions);
      
      rich[t], sprich[t], turnover[t], res_overlap[t], conn[t] = dynstructure(cid,cid_old,sp_v,a_b,tind_m);      

    end
    
    return(cid,rich,sprich,turnover,res_overlap,conn,prim_ext,sec_ext,status,lpot_col)
    
end