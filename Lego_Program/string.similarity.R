string.similarity <- function(x,y) {
  st1 <- x
  st2 <- y
  
  l.st1 <- length(st1)
  l.st2 <- length(st2)
  
  str.intersect <- length(intersect(st1,st2))
  
  Jaccard <- str.intersect / length(unique(c(st1,st2)))
    
  return(Jaccard)
  
}