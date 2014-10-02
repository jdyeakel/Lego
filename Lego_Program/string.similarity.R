string.similarity <- function(x,xab,y,yab) {
  st1 <- x
  st2 <- y
  
  st1.up <- c(st1,setdiff(st2,st1))
  st2.up <- c(st2,setdiff(st1,st2))
  
  ab1 <- xab
  ab2 <- yab
  
  ab1 <- c(ab1,rep(0,length(setdiff(st2,st1))))
  ab2 <- c(ab2,rep(0,length(setdiff(st1,st2))))
  
  #Need to include the organization and weighting of abundances!
  
  m1.sort <- m1[order]
  
  l.st1 <- length(st1)
  l.st2 <- length(st2)
  
  str.intersect <- length(intersect(st1,st2))
  
  Jaccard <- str.intersect / length(unique(c(st1,st2)))
    
  return(Jaccard)
  
}