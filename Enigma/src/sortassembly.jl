function sortassembly(measure,bins,seq)
    reps = size(measure)[1];
    lfseq = findall(x->x>1,diff(seq))[1];
    #This is the initial assembly process
    init_measure = measure[:,1:lfseq];
    
    #this will use reps with >20 nonzero entries
    keepreps = findall(x->x>20,vec(sum(init_measure .> 0,dims=2)));
    
    init_measure_trim = zeros(Float64,length(keepreps),lfseq);
    
    global maxstartmeasure = 0;
    for r=1:length(keepreps)
        # print(r)
        nonzeromeasures = findall(!iszero,init_measure[keepreps[r],:]);
        if length(nonzeromeasures) == 0
            startmeasure = 1;
        else
            startmeasure = nonzeromeasures[1];
        end
        endmeasure = length(init_measure[keepreps[r],:]);
        # init_measure_rm = init_measure[r,findall(!iszero,init_measure[r,:])];
        init_measure_trim[r,1:(endmeasure-startmeasure + 1)] = init_measure[keepreps[r],startmeasure:length(init_measure[keepreps[r],:])]; 
        # init_measure[r,startmeasure:endmeasure];
        #What is the largest number of spots that are skipped?
        global maxstartmeasure = maximum([startmeasure,maxstartmeasure]);
    end
    
    lfseqtrim = maxstartmeasure;
    
    initsteps = bins[bins.<lfseqtrim]; #use these locations for init
    laststeps = bins[bins.>=lfseq]; #use these locations for the rest
    lastbins = indexin(laststeps,seq);

    #Stitch together
    seq_stitch = [initsteps;laststeps];
    measure_stitch = [init_measure_trim[:,initsteps] measure[keepreps,lastbins]];
    
    return measure_stitch, seq_stitch
    
end
