pca_gca <- function(DATA, Jk, cor_min){
  #see Mage, Naes, Hankemeier, & Bro, 2016

  if(missing(cor_min)){
    cor_min <- .9
  }

  num_block <- length(Jk)

  data_block <- list()
  svd_block <- list()
  num_componentBlock <- array() #number of components to be kept per block
  compScores_block <- list()
  compScores_columnlist <- list()

  L <- 1
  for (k in 1:num_block){

    U <- L + Jk[k] - 1
    data_block[[k]] <- DATA[, L:U]
    svdblock <- svd(data_block[[k]])
    cat(sprintf("The eigenvalues of block \"%s\" is", k))
    print(svdblock$d)
    y <- readline("I would like to see the scree plot for the eigenvalues. 1: yes; 0: no. (Please enter 1 or 0.)   ")
    if (y == 1){
      plot(as.vector(svdblock$d), type='b', ylab = "Eigenvalue", xlab = 'Component Number')
    } else if (y!=1 & y!=0){
      stop("Please enter 1 or 0!")
    }
    x <- readline("How many components to be kept for this block?    ")
    num_componentBlock[k] <- as.numeric(x)
    if(num_componentBlock[k]%%1!=0 | num_componentBlock[k]<=0){
      stop("The number of components to be kept must be a positibe integrers!")
    }
    svd_block[[k]] <- svd(data_block[[k]], num_componentBlock[k], num_componentBlock[k])
    compScores_block[[k]] <- svd_block[[k]]$u
    compScores_columnlist <- c(compScores_columnlist, split(compScores_block[[k]], rep(1:ncol(compScores_block[[1]]))))

    L <- U + 1

  }

  #----cononical correlation via rgcca
  compScores_blockCopy <- compScores_block
  compScores_columnlistCopy <- compScores_columnlist
  combinations <- combn(num_block,2)  #this is to generate the index for calculating correlations

  cc <- matrix(1, length(compScores_columnlistCopy), length(compScores_columnlistCopy))
  diag(cc) <- 0
  result.rgcca <- RGCCA::rgcca(compScores_columnlistCopy, cc, tau=rep(0, length(compScores_columnlistCopy)), verbose = FALSE)

  YcompScores_block <- list()
  L <- 1
  for (k in 1:num_block){

    U <- L + num_componentBlock[k] - 1
    YcompScores_block[[k]] <- matrix(unlist(result.rgcca$Y[L:U]), ncol=num_componentBlock[k])
    L <- U + 1
  }

  comdist_indicator <- matrix(0, num_block, max(num_componentBlock))
  complement_template <- matrix(0, num_block, num_block)
  diag(complement_template) <- 1
  complement_indicator <- list()
  remaining_numcomponentBlock <- num_componentBlock

  for (coln in 1: max(num_componentBlock)){

    removeBlock <- array(0) #this is an index for identifitying the blocks that to be removed because its conlumn = 0, and therefore cant enter rgcca.

    for (i in 1: ncol(combinations)){
      if (i %in% removeBlock){
        next
      } else (cor(YcompScores_block[[combinations[,i][1]]][, 1], YcompScores_block[[combinations[,i][2]]][, 1]) >= cor_min){
        comdist_indicator[combinations[,i][1], coln] <- 1
        comdist_indicator[combinations[,i][2], coln] <- 1
      }
   }

    complement_indicator[[coln]] <- complement_template[, which(comdist_indicator[, coln]==0)] #this is to add an extra column for the distinctive component

    if (sum(comdist_indicator[, coln])==0){
      # this means that the entire column is distincive, then the entire procedure can stop.

      for (j in (coln+1):max(num_componentBlock)){
        complement_indicator[[j]] <- complement_template[, which(comdist_indicator[, j]==0)] #this is to add an extra column for the distinctive component
      }

      break

    } else{
      # this means that for column coln, at least two blocks are common. Thus, the procedure continues, and the first column will be removed.
        compScores_columnlist_new <- list()
        for(k in 1:num_block){

          compScores_blockCopy[[k]] <-matrix(compScores_blockCopy[[k]][, -1],nrow = nrow(compScores_blockCopy[[k]]))  #when only a vector left, it must be matrix of Nx1

          if(sum(compScores_blockCopy[[k]])==0) {
            removeBlock <- c(removeBlock, k)
          } else{
            compScores_columnlist_new <- c(compScores_columnlist_new, split(compScores_blockCopy[[k]], rep(1:ncol(compScores_blockCopy[[1]]))))
          }
        }
        remaining_numcomponentBlock <- remaining_numcomponentBlock - 1 #in line with the above for loop, the first column is removed.

        #--- rgcca for the remaining columns
        cc <- matrix(1, length(compScores_columnlist_new), length(compScores_columnlist_new))
        diag(cc) <- 0
        result.rgcca <- RGCCA::rgcca(compScores_columnlist_new, cc, tau=rep(0, length(compScores_columnlist_new)), verbose = FALSE)

        YcompScores_block <- list()
        L <- 1
        for (k in 1:num_block){
          if (remaining_numcomponentBlock[k]>=1){
            U <- L + remaining_numcomponentBlock[k] - 1
            YcompScores_block[[k]] <- matrix(unlist(result.rgcca$Y[L:U]), ncol=num_componentBlock[k])
            L <- U + 1
          } else
            YcompScores_block[[k]] <- 0  #this is to preserve the block structure.
        }
    }

  }


  if(length(complement_indicator)!=0){

    compliment_new <- matrix(NA, nrow=nrow(complement_indicator[[1]]))
    for(m in 1:length(complement_indicator)){
      compliment_new <- cbind(compliment_new, complement_indicator[[m]])
    }
    compliment_new  <- compliment_new [,-1] # remove the NA column
    comdinst_final <- cbind(comdist_indicator, compliment_new)
  }

  comdinst_final <- comdinst_final[, -which(colSums(comdinst_final)==0)]
  for(i in nrow(comdinst_final)){
     comdinst_final <- comdinst_final[, -which(cumsum(comdinst_final[i, ]) > num_componentBlock[i, ])]
  }

  return(comdinst_final)
}
