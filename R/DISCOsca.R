#################################################################################################
# DISCOSCA
# The original function is REVScript_DISCOSCA_TIVvsLAIV.m (obtained from Katrijn Van Deun)
# adjusted by Zhengguo Gu
#################################################################################################

DISCOsca <- function(DATA, R, Jk){



  num_block <- length(Jk)

  pre_results <- svd(DATA, R, R)
  Tmat <- pre_results$u
  Pmat <- pre_results$v %*% diag(pre_results$d[1:R])


  L <- 1
  VAF_results_component <- matrix(NA, num_block, R)
  VAF_results_block <- array(NA, length(Jk))
  ssq_block <- array(NA, length(Jk))

  for (k in 1:length(Jk)){

    U <- L + Jk[k] - 1
    Pmat_k <- Pmat[L:U, ]
    DATA_k <- DATA[, L:U]
    X_hat <- Tmat%*%t(Pmat[L:U, ])
    VAF_results_block[k] <- 1 - sum((X_hat - DATA_k)^2)/sum(DATA_k^2)
    ssq_block[k] <- sum(colSums(Pmat_k^2))
    for (r in 1: R){
      x_pwise_hat <- Tmat[, r] %*% t(Pmat_k[, r])
      VAF_results_component[k, r] <- sum(x_pwise_hat^2)/sum(DATA_k^2)
    }
    L <- U + 1

  }


#find all the combinations

  propExp_Rcomponent <- list() # to record component-wise VAF per block after rotation
  propExp_Rblock <- list() # to record VAF per block after rotation

  TARGET_list <- list()
  TROT_list <- list()
  PROT_list <- list()
  Posit_indicatorList <- list()
  ssq_r_block <- list()

  for(i in max(num_block, R): (R*num_block)){

    combinations <- combn(matrix(1:(R*num_block), num_block, R), i)

    for(j in 1:dim(combinations)[2]){

      posit_indicator <- matrix(0, num_block, R)
      posit_indicator[combinations[,j]] <- 1

      if (min(rowSums(posit_indicator))!=0 & min(colSums(posit_indicator))!=0){
        # if function is to remove cases where an entire column/block is 0

        Posit_indicatorList <- c(Posit_indicatorList, list(posit_indicator))
        W <- matrix(0, sum(Jk), R)

        L <- 1
        for (k in 1:length(Jk)){

          U <- L + Jk[k] - 1
          W_part <- W[L:U, ]
          W_part[, (posit_indicator[k, ]==1)] <- 1
          W[L:U, ] <- W_part
          L <- U + 1

        }

        Target <- 1 - W
        TARGET_list <- c(TARGET_list,list(Target))

        maxiter <- 600;
        convergence <- 0.0000001;
        nrstarts <- 5  #pre=5

        LOSS <- array()
        BMAT <- list()
        for (i in 1:nrstarts){

          B <- pstr(Pmat,Target,W,maxiter,convergence)
          loss <- pstrLoss(B$Bmatrix,Pmat,Target,W)
          LOSS[i] <- loss
          BMAT[[i]] <- B$Bmatrix

        }
        k <- which(LOSS == min(LOSS))
        if (length(k)>1){  #ask Katrijn
          k <- 1
        }
        B <- BMAT[[k]]
        Trot <- Tmat %*% B
        Prot <- Pmat %*% B
        B <- diag(R)

        if(sum(duplicated(posit_indicator,MARGIN=2)) > 0){
            #if > 0, some of the components are the same (in terms of the comm/dist components in the blocks)

            ind <- c(1:R)

            posit_indicator_dup <- posit_indicator

            while(length(ind)>1){

              ind_duplicate <- vector()

              for(j in 2:length(ind)){

                if(sum(posit_indicator_dup[,1] == posit_indicator_dup[, j])==num_block){
                  ind_duplicate <- cbind(ind_duplicate, j)
                }

              }

              if(length(ind_duplicate)==0){
                # in this case, it means either the first column (indexed by ind) is unique, or is the only column left
                ind <- ind[-1]
                posit_indicator_dup <- as.matrix(posit_indicator_dup[, -1]) #as.matrix is needed for the case where only one column is left.
              } else{

              ind_duplicate <- c(1, ind_duplicate) #now we have an index of duplicate columns that are the same as the first column
              ind_dupMATCHind <- ind[ind_duplicate] # here we know which columns in Pmat are duplicated
              dupli_rotate <- svd(Prot[, ind_dupMATCHind], length(ind_dupMATCHind),length(ind_dupMATCHind) )
              B[ind_dupMATCHind, ind_dupMATCHind] <- dupli_rotate$v  #ask Katrijn
              ind <- ind[-ind_duplicate] # remove the duplicated and start searching duplicates in the next columns
              posit_indicator_dup <- as.matrix(posit_indicator_dup[, -ind_duplicate]) #as.matrix is needed for the case where only one column is left.
              }

            }

          }

        Trot <- Trot %*% B
        Prot <- Prot %*% B

        TROT_list <- c(TROT_list, list(Trot))
        PROT_list <- c(PROT_list, list(Prot))

        L <- 1
        ssq_rblock <- array(NA, length(Jk))
        for (k in 1:length(Jk)){

          U <- L + Jk[k] - 1
          Prot_block <- Prot[L:U, (which(posit_indicator[k, ]==1))]
          ssq_rblock[k] <- sum(Prot_block^2)

          L <- U + 1

        }

        ssq_r_block <- c(ssq_r_block, list(ssq_rblock))

        L <- 1
        VAF_rot_component <- matrix(NA, num_block, R)
        VAF_rot_block <- array(NA, length(Jk))
        for (k in 1:length(Jk)){

          U <- L + Jk[k] - 1
          Prot_k <- Prot[L:U, ]

          X_rhat <- Trot%*%t(Prot_k)
          DATA_k <- DATA[, L:U]

          VAF_rot_block[k] <- 1 - sum((X_rhat - DATA_k)^2)/sum(DATA_k^2)

          for (r in 1: R){
            xrot_pwise_hat  <- Trot[, r] %*% t(Prot_k[, r])
            VAF_rot_component[k, r] <- sum(xrot_pwise_hat^2)/sum(DATA_k^2)
          }
          L <- U + 1
        }

        propExp_Rblock <- c(propExp_Rblock, list(VAF_rot_block))
        propExp_Rcomponent <- c(propExp_Rcomponent, list(VAF_rot_component))
      }

    }

  }

  distance <- array(NA, length(ssq_r_block))

  for (j in 1:length(ssq_r_block)){

    distance[j] <- sum((ssq_r_block[[j]]-ssq_block)^2)

  }

  k <- which(distance == min(distance))
  Trot_best <- list(TROT_list[[k]])
  Prot_best <- list(PROT_list[[k]])


  results <- list()

  results$Trot_best <- Trot_best
  results$Prot_best <- Prot_best
  results$k <- c(list(k), min(distance), list(distance), list(ssq_r_block))
  results$propExp_pre_component <- VAF_results_component
  results$propExp_pre_block <- VAF_results_block
  results$propExp_Rotblock <- propExp_Rblock
  results$propExp_Rotcomponent <- propExp_Rcomponent
  results$Target_matrix <- TARGET_list
  results$Postion_indicator <- Posit_indicatorList
  results$Tmatrix <- Tmat
  results$Pmatrix <- Pmat
  results$Trot <- TROT_list
  results$Prot <- PROT_list
  return(results)

}
