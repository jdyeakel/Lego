function sortassembly(measure,bins,seq)
    reps = size(measure)[1];
    lfseq = findall(x->x>1,diff(seq))[1];
    #This is the initial assembly process
    init_measure = measure[:,1:lfseq];
    init_measure_trim = zeros(Float64,reps,lfseq);
    for r=1:reps
        init_measure_rm = init_measure[r,findall(!iszero,init_measure[r,:])];
        init_measure_trim[r,1:length(init_measure_rm)] = init_measure_rm;
    end
    
    initsteps = bins[bins.<lfseq]; #use these locations for init
    laststeps = bins[bins.>=lfseq]; #use these locations for the rest
    lastbins = indexin(laststeps,seq);

    #Stitch together
    seq_stitch = [initsteps;laststeps];
    measure_stitch = [init_measure_trim[:,initsteps] measure[:,lastbins]];
    
    return measure_stitch, seq_stitch
    
end
