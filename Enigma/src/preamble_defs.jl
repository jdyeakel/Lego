function preamble_defs(int_m)

    #Boolian matrices
    a_b = (int_m .== 'a')*1;
    n_b = (int_m .== 'n')*1;
    i_b = (int_m .== 'i')*1;
    m_b = (int_m .== 'm')*1;

    #copy of need binary matrix with diag = 0
    n_b0 = copy(n_b);
    n_b0[diagind(n_b0)] .= 0;

    #Vector and length of species IDs (over all int_m)
    sp_v = findall(isodd,diag(n_b));
    l_sp = length(sp_v);

    int_id = collect(1:size(int_m)[1]);

    return(
    a_b,
    n_b,
    i_b,
    m_b,
    n_b0,
    sp_v,
    int_id
    )
end
