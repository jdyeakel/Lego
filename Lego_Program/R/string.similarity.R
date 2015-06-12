string.similarity <- function(x,xab,y,yab) {
  st1 <- x
  st2 <- y
  
  st1.full <- c(st1,setdiff(st2,st1))
  st2.full <- c(st2,setdiff(st1,st2))
  
  ab1 <- xab
  ab2 <- yab
  
  ab1.full <- c(ab1,rep(0,length(setdiff(st2,st1))))
  ab2.full <- c(ab2,rep(0,length(setdiff(st1,st2))))
  
  #Now the lists have the same elements. Arrange, order, and combine
  st1.ord <- st1.full[order(st1.full)]
  ab1.ord <- ab1.full[order(st1.full)]
  
  st2.ord <- st2.full[order(st2.full)]
  ab2.ord <- ab2.full[order(st2.full)]
  
  if (all(st1.ord == st2.ord) == FALSE) {
    stop("The values are not ordered correctly!")
  } else {
    
    #Need to include the organization and weighting of abundances!
    #ab1.ord <- ab1.ord / sum(ab1.ord)
    #ab2.ord <- ab2.ord / sum(ab2.ord)
    
    #db <- data.frame(st1.ord,ab1.ord,ab2.ord)
    
    #Cosine similarity index
    cosine.sim <- (ab1.ord %*% ab2.ord) / sqrt(((ab1.ord %*% ab1.ord) * (ab2.ord %*% ab2.ord)))
    
    str.intersect <- length(intersect(st1,st2))
    
    Jaccard <- str.intersect / length(unique(c(st1,st2)))
    
    output <- c(Jaccard,cosine.sim)
    
    return(output)
    
  }
  
}