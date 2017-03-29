#------------- bound function --------------
.bound <- function(x, bds){
  x[x > max(bds)] <- max(bds)
  x[x < min(bds)] <- min(bds)
  x
}
