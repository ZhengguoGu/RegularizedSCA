#################################################################################################
# DISCOSCA
# The original function is REVScript_DISCOSCA_TIVvsLAIV.m (obtained from Katrijn Van Deun)
# adjusted by Zhengguo Gu
#################################################################################################

DISCOsca <- function(Pmat, R, DATA, Jk){

  num_block <- length(Jk)

#find all the combinations

  for(i in 2: (R*num_block)){

    combinations <- combn(matrix(1:(R*num_block), num_block, R), i)

    for(j in 1:dim(combinations)[2]){

      posit_indicator <- matrix(0, num_block, R)
      posit_indicator[combinations[,j]] <- 1

      if (min(rowSums(posit_indicator))==0 | min(colSums(posit_indicator))==0){

        break

        } else {

          W <- matrix(0, sum(Jk), R)
          L <- 1
          for (k in 1:length(Jk)){

            U <- L + Jk[i] - 1
            P_postion <- W[L:U, ]
            P_position[, posit_indicator[k, ]] <- 1
            L <- U + 1

          }

          Target <- 1 - W

          maxiter <- 600;
          convergence <- 0.0000001;
          nrstarts <- 5

          LOSS <- array()
          BMAT <- list()
          for (i in 1:nrstarts){

            B <- pstr(Pmat,Target,W,maxiter,convergence)
            loss <- pstrLoss(B$Bmatrix,Pmat,Target,W)
            LOSS[i] <- loss
            BMAT[[i]] <- B$Bmatrix

          }
          k <- which(LOSS == min(LOSS))
          B <- BMAT[[k]]
        }

    }

  }
  zz <- combn(dd, 4)
}
