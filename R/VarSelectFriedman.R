#variable selection with lasso and group lasso. The algorithm is developed based on Friedman, Hastie, and Tibshirani (2010)

VarSelectFriedman <- function(DATA, Jk, R, LASSO, GROUPLASSO, MaxIter){

  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)

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
        r <- as.vector(DATA[ ,c(L:U)]) # Note here, r is a vector and is different from Friedman et al.

        sum_abs_theta <- sum(abs(Pt_1))
        if (sum_abs_theta == 0){
           break  #this is because the P's are already zeros
        }
        else{

          Zij <- kronecker(diag(Jk[i]), Tmat)
          X_k_r <- t(Zij) %*% r
          soft_Xkrk_lumba <- array()
          for (k in 1: length(X_k_r)){

            soft_Xkrk_lumba[k] <- sign(X_k_r[k])*max((abs(X_k_r[k]) - LASSO), 0)

          }

          if (sqrt(sum(soft_Xkrk_lumba^2)) / GROUPLASSO^2 <= 1) {

            Pt[ ,c(L:U)] <- 0
            P <- t(Pt)

          }

          else {

            Pt_temporaryVec <- as.vector(Pt_1)

            for (j in 1:length(Pt_temporaryVec)){

              Pt_temporaryVec[j] <- 0 # remove the variable that is to be estimated
              Pt_temporary <- matrix(Pt_temporaryVec, R, Jk[i])
              wj <- r - as.vector(Zij %*% Pt_temporaryVec)

              if((t(Zij[, j]) %*% wj) < LASSO){

                Pt[ ,c(L:U)] <- Pt_temporary
                P <- t(Pt)
              }
              else{
                f <- function(x) {
                  0.5*sum((wj - Zij[,j]*x)^2) + GROUPLASSO*sqrt(Jk[i])*(sum((Pt_temporary^2))+x^2)^0.5 + LASSO*(sum(abs(Pt_temporary))+abs(x))
                  }
                xmin <- optimize(f, c(-1,1), tol = 0.0001)
                Pt_temporaryVec[j] <- xmin$minimum
                Pt[ ,c(L:U)] <- matrix(Pt_temporaryVec, R, Jk[i])
                P <- t(Pt)
              }
            }
          }
        }
        L <- U + 1
      }
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
  return_varselect$residual <- residual
  return_varselect$lassopen <- lassopen
  return_varselect$glassopen <- Glassopen
  return_varselect$iter <- iter - 1
  return(return_varselect)


}
