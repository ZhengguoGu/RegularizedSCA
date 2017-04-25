
#Variable selection with Lasso and Group Lasso - Version 1: "datablock" method. 
#The algorithm applies sparse Group Lasso on the entire concatenate datablock, 

CDfriedmanV1 <- function(DATA, Jk, R, LASSO, GROUPLASSO, MaxIter){
  
  DATA <- data.matrix(DATA)
  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)

  if(missing(MaxIter)){
    MaxIter <- 400
  }

  #initialize P
  P <- matrix(rnorm(sumJk * R), nrow = sumJk, ncol = R)
  Pt <- t(P)

  absP <- abs(P)
  pen_l <- LASSO * sum(absP)
  sqP <- P^2
  L <- 1
  pen_g <- 0
  for (i in 1:length(Jk)){

    U <- L + Jk[i] - 1
    sqrtsumP <- sqrt(colSums(sqP[L:U, ]))/sqrt(Jk[i])
    pen_g <- pen_g + GROUPLASSO * sum(sqrtsumP) * Jk[i]
    L <- U + 1
  }

  residual <- sum(DATA ^ 2)
  Lossc <- residual + pen_l + pen_g

  conv <- 0
  iter <- 1
  Lossvec <- array()

  while (conv == 0){

    #update Tmat, note that Tmax refers to T matrix
    if (LASSO == 0 & GROUPLASSO == 0){
      SVD_DATA <- svd(DATA, R, R)  #note this is different from the matlab svds function. need to test it!!
      Tmat <- SVD_DATA$u
    }
    else {
      A <- Pt %*% t(DATA)
      SVD_DATA <- svd(A, R, R)
      Tmat <- SVD_DATA$v %*% t(SVD_DATA$u)
    }

    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu <- residual + pen_l + pen_g

    #update P
    if (LASSO == 0 & GROUPLASSO == 0){

      P <- t(DATA) %*% Tmat
      Pt <- t(P)

    }
    else{

      L <- 1
      for (i in 1:length(Jk)){ #iterate over groups

        U <- L + Jk[i] - 1
        Pt_1 <- Pt[ ,c(L:U)]
        data <- DATA[ ,c(L:U)] #here we focus on each datablock.

        if (sum(abs(Pt_1)) != 0){
         
          #  calculate S(2Vec(T'X_l), lambda_L), and its l2 norm.
          
          S_2Vec_Lambda <- soft_th(2 * t(Tmat) %*% data, LASSO)
          
          S_2Vec_Lambda_norm <- sqrt(sum(S_2Vec_Lambda^2))
          
          s_l2 <- 0.5 - 0.5 * GROUPLASSO * sqrt(Jk[i])/S_2Vec_Lambda_norm
         
          if (s_l2 <= 0){
        
            Pt[ ,c(L:U)] <- 0
          
          } else{
            Pt[ ,c(L:U)] <- s_l2 * S_2Vec_Lambda
          }
         }
        
        L <- U + 1
      }
      P <- t(Pt)
    }


    pen_l <- LASSO*sum(abs(P))
    sqP <- P^2
    L <- 1
    pen_g <- 0
    for (i in 1:length(Jk)){
      U <- L + Jk[i] - 1
      sqrtsumP <- sqrt(colSums(sqP[L:U, ]))/sqrt(Jk[i])
      pen_g <- pen_g + GROUPLASSO * sum(sqrtsumP) * Jk[i]
      L <- U + 1
    }

    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu2 <- residual + pen_l + pen_g

    #check convergence
    if (abs(Lossc-Lossu)< 10^(-9)) {
      Loss <- Lossu
      residual <- residual
      lassopen <- pen_l
      Glassopen <- pen_g
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    else if (iter > MaxIter | LASSO == 0){
      Loss <- Lossu
      residual <- residual
      lassopen <- pen_l
      Glassopen <- pen_g
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
  #return_varselect$glassopen <- Glassopen
  #return_varselect$iter <- iter - 1
  return(return_varselect)

}
