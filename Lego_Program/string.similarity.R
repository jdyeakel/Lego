string.similarity <- function(x,xab,y,yab) {
  st1 <- x
  st2 <- y
  
  st1.ab <- xab
  st2.ab <- yab
  
  l.st1 <- length(st1)
  l.st2 <- length(st2)
  
  is.same <- numeric(l.st1)
  for (i in 1:l.st1) {
    if (st1[i] %in% st2) {
      is.same[i] <- 1
    } else {is.same[i] <- 0}
  }
  tot.same <- sum(is.same)
  
  overlap <- which(is.same==1)
  
  ab.diff <- 
  
}