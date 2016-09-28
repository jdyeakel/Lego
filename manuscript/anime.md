# Description 

*	Iteratively build a template - how many different communities emerge?
* What does this depend on? (number of primary producers, generality of founding species, etc
*	Dependency of `n` interaction vs. template size
*	Sensitivity of *anime* communities vs. traditional trophic networks (with `n` and `m` at zero probability)


### Building full communities over template size, and threshold variations  
*	Template size, threshold value  
*	Determine richness, variation in community assembly process (by what metric???) 
*	Average this over multiple repetitions for various `int_m`  
*	

### Colonization Module  
*	Should it be random? For now, yes  
*	Integrate competition check (based on similarity)  
*	Similarity: `sum(seq1 .== seq2)/num_play`  This will be a `num_play x num_play` matrix with the diagonal = 1 to check  
*	

### Extinction Module  
*	Randomly select species, or construct an extinction risk vector that is based on: predation load, competition (similarity)  
*	


Invasive species idea: in spatial community simulation, after regions have formed, transplant a species from one habitat to another to determine impact  
