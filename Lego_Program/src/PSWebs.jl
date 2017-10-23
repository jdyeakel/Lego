function PSWebs(com_tind,com_mind,sigma,reps)
    
    Jm = Jmatrix(com_tind,com_mind);
    AJm = abs(Jm);
    # R"image($Jm,col=gray.colors(2))"
    
    #stability analysis over sigma
    maxre = zeros(reps);
    @sync @parallel for r=1:reps
        #create random matrix
        intdist = Normal(0,sigma);
        Rm = rand(intdist,size(Jm,1),size(Jm,2));
        Jacobian = Rm .* Jm;
        #reset the diagonal to -1
        for i=1:size(Jm,1)
            Jacobian[i,i] = -1;
        end
        evalue = eig(Jacobian);
        re = real.(evalue)[1];
        ie = imag.(evalue)[1];
        maxre[r] = maximum(re);
        # R"plot($re,$ie)"
    end
    PSW = length(find(x->x<1*10^-6,maxre))/reps;

    return(PSW);

end