# Variable selection algorithm with a predefined component loading structure

CDpre <- function(DATA, Jk, R, CommPosition, GroupStructure, LASSO, MaxIter){

  DATA <- data.matrix(DATA)
  DistPosition <- setdiff(1:R, CommPosition)
  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)

  if(missing(MaxIter)){
    MaxIter <- 400
  }

  #initialize P, Lossc
  P <- matrix(rnorm(sumJk * R), nrow = sumJk, ncol = R)
  P[GroupStructure == 0]<-0
  Pt <- t(P)

  PIndexforLasso <- Pt
  PIndexforLasso[CommPosition, ] <- 1
  PIndexforLasso[DistPosition, ] <- 0
  PIndexforGLasso <- Pt #note that the distinctive structure is predefined.
  PIndexforGLasso[CommPosition, ] <- 0
  PIndexforGLasso[DistPosition, ] <- 1

  pen1 <- LASSO*sum(abs(P[, CommPosition]))
  sqP <- P^2
  L <- 1

  residual <- sum(DATA^2)
  Lossc <- residual + pen1

  conv <- 0
  iter <- 1
  Lossvec <- array()

  while (conv == 0){

    #update Tmat, note that Tmax refers to T matrix
    if (LASSO == 0){
      SVD_DATA <- svd(DATA, R, R)
      Tmat <- SVD_DATA$u
    }
    else {
      A <- Pt %*% t(DATA)
      SVD_DATA <- svd(A, R, R)
      Tmat <- SVD_DATA$v %*% t(SVD_DATA$u)
    }

    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu <- residual + pen1

    #update P
    if (LASSO == 0){
      P <- t(DATA) %*% Tmat
      P[GroupStructure == 0] <- 0  # this is to keep the zero structures
      Pt <- t(P)
    }
    else{

      for (r in 1:R){
        if (r %in% CommPosition) {
          for (j in 1:sumJk){
            ols <- t(DATA[, j]) %*% Tmat[, r]
            Lambda <- 0.5 * LASSO
            if (ols < 0 & abs(ols) > Lambda) {
              P[j, r] <- ols + Lambda
            }
            else if (ols > 0 & abs(ols) > Lambda) {
              P[j, r] <- ols - Lambda
            }
            else {
              P[j, r] <- 0
            }
          }
        }
        else {
          for (j in 1:sumJk){
            P[j, r] <- t(DATA[, j]) %*% Tmat[, r] #note that in the original matlab file this term is devided by sumD(=1)
          }
        }
      }
    }

    P[GroupStructure == 0] <- 0
    Pt <- t(P)

    #absP <- abs(P)
    pen1 <- LASSO*sum(abs(P[, CommPosition]))
    sqP <- P^2

    L <- 1

    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu2 <- residual + pen1

    #check convergence
    if (abs(Lossc-Lossu)< 10^(-9)) {
      Loss <- Lossu
      residual <- residual
      lassopen <- pen1
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    else if (iter > MaxIter | LASSO == 0){
      Loss <- Lossu
      residual <- residual
      lassopen <- pen1
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }

    Lossvec[iter] <- Lossu
    iter <- iter + 1
    Lossc <- Lossu2
  }

  return_varselect <- list()
  return_varselect$Pmatrix <- P
  return_varselect$Tmatrix <- Tmat
  return_varselect$Loss <- Loss
  return_varselect$Lossvec <- Lossvec
  #return_varselect$residual <- residual
  #return_varselect$lassopen <- lassopen
  #return_varselect$iter <- iter - 1
  return(return_varselect)

}
