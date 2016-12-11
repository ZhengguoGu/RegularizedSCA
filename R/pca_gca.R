pca_gca <- function(DATA, Jk, eig_min, cor_min){
  #see Mage, Naes, Hankemeier, & Bro, 2016
  num_block <- length(Jk)

  data_block <- list()
  svd_block <- list()
  num_componentBlock <- array() #number of components to be kept per block
  L <- 1
  for (k in 1:num_block){

    U <- L + Jk[k] - 1
    data_block[[k]] <- DATA[, L:U]
    svd_block[[k]] <- svd(data_block[[k]])
    cat(sprintf("The eigenvalues of block \"%s\" is", k))
    print(svd_block[[k]]$d)
    print("And the screeplot is")
    plot(as.vector(svd_block[[k]]$d), type='b', ylab = "Eigenvalue", xlab = 'Component Number')
    x <- readline(promt = "How many components to be kept for this block? (Please fill in a positive integer, such as 1, 2, 3, etc.)")
    num_componentBlock[k] <- as.numeric(x)
    if(is.na(x) | x<=0){
      stop("The number of components to be kept must be a positibe integrers!")
    }


    L <- U + 1

  }





}
