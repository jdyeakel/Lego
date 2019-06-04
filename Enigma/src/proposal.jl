function proposal(temperature,value,valuerange)
    newrange_min = Int64(maximum([findall(x->x==value,valuerange)[1]-Int64(floor((length(valuerange)*temperature)/2)),1]));
    newrange_max = Int64(minimum([findall(x->x==value,valuerange)[1]+Int64(floor((length(valuerange)*temperature)/2)),length(valuerange)]));
    newrange = valuerange[newrange_min:newrange_max];
    
    #select the new proposal
    
    proposal_value = rand(newrange);
    
    return proposal_value
end