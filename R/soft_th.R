## function soft-thresholding

soft_th <- function(X, lambda){
  
  # assume X is a matrix
  result <- X
  index1 <- which(X > lambda)
  result[index1] <- X[index1] - lambda
  index2 <- which(X < -lambda)
  result[index2] <- X[index2] + lambda
  index0 <- which(X <= lambda & X >= -lambda)
  result[index0] <- 0
  
  return(result)
}
