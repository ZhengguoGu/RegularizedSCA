#########################################################################
# pstr.R finds a (orthogonal) rotation to a partially specified target
#
# The following objective function is used:
#      min_B ||Wo(PB-Target)||2
# with W a binary matrix of weights
#      P the matrix to rotate
#      Target the target to reach
#
# Input P: matrix to rotate
#       Target: target to reach
#       W: weight matrix (0/1 for resp. unspecifed and specified elements of the target)
#       maxiter: maximum number of iterations
#       convergence: minimal loss between current and previous iteration
# Output B: an orthogonal rotation matrix
#
# Author: Katrijn Van Deun; Zhengguo Gu implemented in R
##########################################################################

### Ask Katrijn: The B matrices are the name???

pstr <- function(P, Target, W, maxiter,convergence){
  n <- dim(P)[1]
  m <- dim(P)[2]

  L <- array()
  Bmat <- list()
  REFL <- reflexmat(m) # requires reflexmat.R

  for(i in 1: dim(REFL)[1]){

    k <-which(REFL[i, ] == -1)
    Binit <- diag(m)
    Binit[, k] <- -1*Binit[, k]

    B1 <- t(P) %*% P  #ask katrijn
    alpha <- max(eigen(B1)$values)
    iter <- 1

    stop <- 0

    Bcurrent <- Binit
    Lossc <- pstrLoss(Binit, P, Target, W)

    while(stop == 0){

      Pw <- W * Target + P %*% Bcurrent - W * (P %*% Bcurrent)
      A <- -2 * t(Pw) %*% P
      Fmat <- A + 2 * t(Bcurrent) %*% t(B1) - 2 * alpha * t(Bcurrent)
      F_svd <- svd(-Fmat)
      B <- F_svd$v %*% t(F_svd$u)

      if (iter == maxiter){
        stop <- 1
      }

      Loss <- pstrLoss(B, P, Target, W)
      Diff <- Lossc - Loss

      if(abs(Diff) < convergence){
        stop <- 1
      }

      iter <- iter + 1
      Lossc <- Loss
      Bcurrent <- B
    }

    L[i] <- Lossc
    Bmat[[i]] <- Bcurrent
  }

  k <- which(L == min(L))
  Loss <- L[k[1]]
  B <- Bmat[[k[1]]]


  results <- list()
  results$Bmatrix <- B
  results$Loss <- Loss
  return(results)
}



