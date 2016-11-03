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
*	NOTE: need to deal with boarder-case scenarios

> **Ideas**
> Assess fixed point richness, similarity, connectance, other food web metrics?, mutualistic net metric as a function of 1) pr(n), pr(m) for a given pr(a) which will result in a specific connectance/size relationship
>	Replicate primary extinction -> number of secondary extinctions from Sole/Roopnarine -> (!!)encoded in the extinction module (**SAVE DATA**)
>	Perturbation idea: Random perturbation vs. targeted perturbation, where the target is defined by a suite of species extinctions that have SIMILAR interactions
>	Invasive species idea: in spatial community simulation, after regions have formed, transplant a species from one habitat to another to determine impact


### Notes/Problems
* How to account for the number of timesteps that the system spends in a certain state and the potential for autocorrelation? Bootstrapping?
* 
