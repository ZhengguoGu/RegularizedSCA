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
    y <- readline("I would like to see the scree plot for the eigenvalues. 1: yes; 0: no. (Please enter 1 or 0.) ")
    if (y == 1){
      plot(as.vector(svdblock$d), type='b', ylab = "Eigenvalue", xlab = 'Component Number')
    } else if (y!=1 & y!=0){
      stop("Please enter 1 or 0!")
    }
    x <- readline("How many components to be kept for this block?")
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
  c <- matrix(1, length(compScores_columnlist), length(compScores_columnlist))
  diag(c) <- 0
  result.rgcca <- RGCCA::rgcca(compScores_columnlist, c, tau=rep(0, length(compScores_columnlist)), verbose = FALSE)

  combinations <- combn(num_block,2)  #this is to generate the index for calculating correlations
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
  compScores_blockCOPY <- compScores_block  # this copy is going to be used when redo rgcca

  for (c in 1: max(num_componentBlock)){

    redoRgcca <- 0 #if the entire column is considered distictive, then the rest columns will be distinctive.
    for (i in 1: ncol(combinations)){
      if (cor(YcompScores_block[[combinations[,i][1]]][, 1], YcompScores_block[[combinations[,i][2]]][, 1]) >=.cor_min){
        comdist_indicator[combinations[,i][1], c] <- 1
        comdist_indicator[combinations[,i][2], c] <- 1
      }
    }

    complement_indicator[[c]] <- complement_template[, which(comdist_indicator[, c]==0)]

    #---- remove the first column and redo rgcca
    if(sum(comdist_indicator[, c])!=0){

      compScores_columnlist_new <- list()
      for(k in 1:num_block){

        if(ncol(compScores_blockCOPY[[k]])==0){
          compScores_columnlist_new <- c(compScores_columnlist_new, split(compScores_blockCOPY[[k]], rep(1:ncol(compScores_blockCOPY[[1]]))))
        } else{
          compScores_blockCOPY[[k]] <- compScores_blockCOPY[[k]][, -1]
          compScores_columnlist_new <- c(compScores_columnlist_new, split(compScores_blockCOPY[[k]], rep(1:ncol(compScores_blockCOPY[[1]]))))
        }
      }


    }

  }
}
