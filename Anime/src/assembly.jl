function assembly(
    int_m,a_b,n_b,i_b,m_b,n_b0,sp_v,int_id,tp_m,tind_m,
    a_thresh,n_thresh,extinctions,tmax)
    
    cid = Array{Int64}(0);
    rich = Array{Int64}(tmax);
    sprich = Array{Int64}(tmax);
    turnover = Array{Float64}(tmax);
    prim_ext = Array{Int64}(tmax);
    sec_ext = Array{Int64}(tmax);
    res_overlap = Array{Float64}(tmax);
    conn = Array{Float64}(tmax);
    conn_ind = Array{Float64}(tmax);
    status = Array{Int64}(tmax);
    lpot_col = Array{Int64}(tmax);
    avgdegree = Array{Float64}(tmax);
    
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
      
      rich[t], sprich[t], turnover[t], res_overlap[t], conn[t], conn_ind[t] = dynstructure(cid,cid_old,sp_v,a_b,tp_m,tind_m);      
      
      # spcid = intersect(sp_v,cid);
      # spcid_ind = indexin(spcid,[1;sp_v]);
      # degrees = deleteat!(vec(sum(tind_m[[1;spcid_ind],[1;spcid_ind]],2)),1);
      
      #These are 'potential degrees'
      degrees = sum(a_b[cid,:],2);
      avgdegree[t] = mean(degrees);

    end
    
    return(cid,rich,sprich,turnover,res_overlap,conn,conn_ind,prim_ext,sec_ext,status,lpot_col,avgdegree)
    
end