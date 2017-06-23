#' DISCO-SCA rotation.
#'
#' A DISCO-SCA procedure for identifying common and distinctive components.
#'
#' @param DATA A matrix, which contains the concatenated data with the same subjects from multiple blocks.
#' Note that each row represents a subject.
#' @param Jk A vector containing number of variables in the concatenated data matrix.
#' @param R Number of components (R>=2).
#'
#' @return
#' \item{Trot_best}{Estimated component score matrix (i.e., T)}
#' \item{Prot_best}{Estimated component loading matrix (i.e., P)}
#' \item{comdist}{A matrix representing common distinctive components. (Rows are data blocks and columns are components.) 0 in the matrix indicating that the corresponding
#' component of that block is estimated to be zeros, and 1 indicates that (at least one component loading in) the corresponding component of that block is not zero.
#' Thus, if a column in the \code{comdist} matrix contains only 1's, then this column is a common component, otherwise distinctive component.}
#' \item{propExp_component}{Proportion of variance per component.}
#' @examples
#' \dontrun{
#' DATA1 <- matrix(rnorm(50), nrow=5)
#' DATA2 <- matrix(rnorm(100), nrow=5) #thus, we assume that DATA1 and DATA2 are 
#'                                     #with respect to the same 5 subjects here.
#' DATA <- cbind(DATA1, DATA2)
#' R <- 5
#' Jk <- c(10, 20) #DATA1 has 10 columns, DATA2 20.
#' DISCOsca(DATA, R, Jk)
#'}
#' @references
#' Schouteden, M., Van Deun, K., Wilderjans, T. F., & Van Mechelen, I. (2014). 
#' Performing DISCO-SCA to search for distinctive and common information in linked data. 
#' Behavior research methods, 46(2), 576-587.
#' @export

DISCOsca <- function(DATA, R, Jk){

  DATA <- data.matrix(DATA)
  num_block <- length(Jk)

  if(R == 1){
    stop("Parameter R = 1 is not allowed.")
  }
  
  if (R >= min(dim(DATA))){
    stop("# of components must be smaller than # of subjects (i.e., # of rows of the concatenated data) and # of variables (i.e., # of columns of the concatenated data")
  }
  pre_results <- svd(DATA, R, R)
  Tmat <- pre_results$u
  Pmat <- pre_results$v %*% diag(pre_results$d[1:R])


  L <- 1
  VAF_results_component <- matrix(NA, num_block, R)
  VAF_results_block <- array(NA, num_block)
  ssq_block <- array(NA, num_block)
  distance <- array(NA)

  for (k in 1:num_block){

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
  #propExp_Rblock <- list() # to record VAF per block after rotation

  TARGET_list <- list()
  TROT_list <- list()
  PROT_list <- list()
  Posit_indicatorList <- list()
  #ssq_r_block <- list()

  for(i in max(num_block, R): ((R*num_block)-1)){

    combinations <- utils::combn(matrix(1:(R*num_block), num_block, R), i)

    for(j in 1:dim(combinations)[2]){

      posit_indicator <- matrix(0, num_block, R)
      posit_indicator[combinations[,j]] <- 1
      
      if (min(rowSums(posit_indicator))!=0 & min(colSums(posit_indicator))!=0){
        # if function is to remove cases where an entire column/block is 0

        cat(sprintf("Now checking the following component structure:\n"))
        print(posit_indicator)
        
        Posit_indicatorList <- c(Posit_indicatorList, list(posit_indicator))
        Target <- matrix(0, sum(Jk), R)

        L <- 1
        for (k in 1:num_block){

          U <- L + Jk[k] - 1
          Target_part <- Target[L:U, ]
          Target_part[, (posit_indicator[k, ]==1)] <- 1
          Target[L:U, ] <- Target_part
          L <- U + 1

        }

        W <- 1 - Target
        TARGET_list <- c(TARGET_list,list(Target))

        maxiter <- 5000;
        convergence <- 0.0001;
        nrstarts <- 2  #default was 5 in fact, but its too slow. 

        LOSS <- array()
        BMAT <- list()
        for (i in 1:nrstarts){

          B <- pstr(Pmat,Target,W,maxiter,convergence)
          loss <- pstrLoss(B$Bmatrix,Pmat,Target,W)
          LOSS[i] <- loss
          BMAT[[i]] <- B$Bmatrix

        }
        k <- which(LOSS == min(LOSS))
        B <- BMAT[[k[1]]]
        Trot <- Tmat %*% B
        Prot <- Pmat %*% B
        B <- diag(R)


        ####################### rotate entire component across blocks ##############################
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
        ##################################################################################################
        Trot <- Trot %*% B
        Prot <- Prot %*% B

        TROT_list <- c(TROT_list, list(Trot))
        PROT_list <- c(PROT_list, list(Prot))


        L <- 1
        ssq_r_component <- matrix(NA, num_block, R)

        for (k in 1:num_block){

          U <- L + Jk[k] - 1
          for (r in 1: R){

            ssq_r_component[k, r] <- sum(Prot[L:U, r]^2)
          }
          L <- U + 1

        }
        propExp_Rcomponent <- c(propExp_Rcomponent, list(ssq_r_component))

        distance_component <- array(NA, R)
        for (r in 1:R){

          if(sum(posit_indicator[, r])==num_block){
            #if(true), we have a common component
            combi <- utils::combn(1:num_block, 2) #when there are more than 2 blocks, combi is useful

            distance_combi <- array(NA, dim(combi)[2])
            for (c in 1:dim(combi)[2]){
              distance_combi[c] <- abs(ssq_r_component[combi[, c][1], r]/ssq_block[combi[, c][1]]-ssq_r_component[combi[, c][2], r]/ssq_block[combi[, c][2]])
            }

            distance_component[r] <- max(distance_combi) # ask Katrijn

          } else{

            zeros <- which(posit_indicator[, r] ==0)
            distance_zero <- array(NA, length(zeros)) # there may be cases where more than 2 blocks have zeros
            for(z in 1: length(zeros)){
              distance_zero[z] <- ssq_r_component[zeros[z], r]/ssq_block[zeros[z]]
            }
            distance_component[r] <- max(distance_zero)
          }
        }

        distance <- c(distance, max(distance_component))
      }

    }

  }

  distance <- distance[-1] # the first one is always NA

  k <- which(distance == min(distance))
  Trot_best <- list(TROT_list[[k]])
  Prot_best <- list(PROT_list[[k]])


  results <- list()

  results$Trot_best <- Trot_best
  results$Prot_best <- Prot_best
  results$comdist <- Posit_indicatorList[k]
  #results$k <- c(list(k), min(distance), list(Posit_indicatorList[k]))
  #results$propExp_pre_component <- VAF_results_component
  #results$propExp_pre_block <- VAF_results_block
  #results$propExp_Rotblock <- propExp_Rblock
  results$propExp_component <- propExp_Rcomponent[k]
  #results$Target_matrix <- TARGET_list
  #results$Postion_indicator <- Posit_indicatorList
  #results$Tmatrix <- Tmat
  #results$Pmatrix <- Pmat
  #results$Trot <- TROT_list
  #results$Prot <- PROT_list
  return(results)

}
