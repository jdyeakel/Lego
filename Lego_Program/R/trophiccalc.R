species <- LETTERS[1:6];
links <- c(0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
L <- matrix(links, nrow = 6, byrow = TRUE, dimnames = list(FROM=species, TO=species))


library(expm)

calcHeight <- function(MAT) {
  ## Find 'leaf nodes' (i.e. species that are only eaten, 
 ## and don't eat any others) 
 leaves <- which(rowSums(L)==0) 
 ## Find the maximum possible chain length (if this is a DAG) 
 maxHeight <- nrow(MAT) - length(leaves) - 1 
 ## Then use it to determine which matrix powers we'll need to calculate. 
 index <- seq_len(maxHeight) 
 paths <- lapply(index, FUN=function(steps) MAT %^% steps) 
 pathSteps <- lapply(index, FUN=function(steps) (1 + steps) * paths[[steps]]) 
 ## Trophic height is expressed relative to leaf nodes 
 paths <- Reduce("+", paths)[-leaves, leaves] 
 pathSteps <- Reduce("+", pathSteps)[-leaves, leaves] 
 rowSums(pathSteps)/rowSums(paths) 
} 

calcHeight(L)
