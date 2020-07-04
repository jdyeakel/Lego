# ENIgMa ecological assembly model

Hi - thanks for your interest in the ENIgMa ecological assembly model. We recently published our paper describing and analyzing the model [here](https://www.nature.com/articles/s41467-020-17164-x). Within this repository is the code that we used for the paper. I'm hoping to soon clean this up a bit and make it more straightforward, but not easy with a 5 year old pulling on my face all of the time.

The basic assembly model can be run by: 
> 1) loading the functions in `src/`
> 2) running lines 7-49 in `sim.jl` 

Here is the list of the output: 
> `sprich`: vector of species richness over time
> `rich`: vector of total richness over time (species + modifiers)
> `clock`: time vector
> `CID`: boolean array with species/modifier IDs as rows and simulation steps as column. TRUE=present
> `events`: a vector that tracks the different events that occur. Colonization = 0; Primary extinction = 1; Secondary extinction = 2; Modifier extinction = 3.

