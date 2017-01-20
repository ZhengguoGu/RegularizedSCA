#' PCA-GCA method
#'
#' Use PCA-GCA method to identify common and distictive components.
#'
#' @param DATA A concatinated data matrix with the same number of rows.
#' @param Jk A vector containing number of variables  in the concatinated data matrix. Please see the example below.
#' @param cor_min The minimum corelation bewtween two components. The default value is .9; thus, it means that if the correlation
#' between the two component is at least .9, then these two components are regarded as forming a single common component.
#' @return It prints out the number of components of each block and the number of common components. It also returns the component scores for each block for further analysis.

#' @examples
#' DATA <- c(DATA1, DATA2)
#' Jk <- c(ncol(DATA1), ncol(DATA1))
#' pca_gca(DATA, Jk, cor_min = .8)
#'
#' @references Tenenhaus, A., & Tenenhaus, M. (2011). Regularized generalized canonical correlation analysis. Psychometrika, 76(2), 257-284
#' @note
#' Please be ware of the interactive input: The function first performs PCA on each data block and then displays the eigenvalues (and a scree plot).
#' Afterwards the function awaits the input from the user - it needs to know how many components need to be retained for that block.
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

    L <- U + 1

  }

  #----cononical correlation via rgcca

  canonical_cor <- RGCCA::rgcca(compScores_block, C=1-diag(length(compScores_block)), tau = rep(0, length(compScores_block)),
                                ncomp = min(num_componentBlock), verbose = FALSE)

  com_comp <- sum(sqrt(canonical_cor$AVE$AVE_inner) >= cor_min)

  cat("The number of components in each block are:", num_componentBlock)
  cat(sprintf("There are in total %s common components in the concatenated data.", com_comp))

  return(compScores_block)

}
