ipbc <- function(x,L) {
  check <- 0
  #Bottom Row
  if (x <= L+2) {new_x <- x + L*(L+2); check <- 1}
  #Top Row
  if (x >= ((L+2)^2 - (L+1))) {new_x <- x - L*(L+2); check <- 1}
  #Left Row
  if (x%%(L+2) == 1) {new_x <- x + L; check <- 1}
  #Right Row
  if (x%%(L+2) == 0) {new_x <- x - L; check <- 1}
  #Bottom Left
  if ((x <= L+2) && (x%%(L+2) == 1)) {new_x <- x + L*(L+2) + L; check <- 1}
  #Bottom Right
  if ((x <= L+2) && (x%%(L+2) == 0)) {new_x <- x + L*(L+2) - L; check <- 1}
  #Top Left
  if ((x >= ((L+2)^2 - (L+1))) && (x%%(L+2) == 1)) {new_x <- x - L*(L+2) + L; check <- 1}
  #Top Right
  if ((x >= ((L+2)^2 - (L+1))) && (x%%(L+2) == 0)) {new_x <- x - L*(L+2) - L; check <- 1}
  
  #Middle
  if(check == 0) {new_x <- x}
  nn <- c(new_x+1,new_x-1,new_x+(L+2),new_x-(L+2))
  return(nn)
}