pca_gca <- function(DATA, Jk, eig_min, cor_min){
  #see Mage, Naes, Hankemeier, & Bro, 2016
  num_block <- length(Jk)

  data_block <- list()
  svd_block <- list()
  num_componentBlock <- array() #number of components to be kept per block
  component_block <- list()
  Tcomponent_block <- list()
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
    component_block[[k]] <- as.matrix(svd_block[[k]]$v) %*% diag(svd_block[[k]]$d[1:num_componentBlock[k]])
    Tcomponent_block[[k]] <- t(component_block[[k]])
    L <- U + 1

  }

  #----cononical correlation via rgcca
  c <- matrix(1, num_block, num_block)
  diag(c) <- 0
  result.rgcca <- RGCCA::rgcca(data_block, c, tau=rep(0, num_block), verbose = FALSE)
  ax <- list()
  for (l in 1:num_block){
    ax[[l]] <- sweep(component_block[[l]], 2, result.rgcca$a[[l]], "*")
  }
  cor(ax[[1]][,3], ax[[2]][,4])
  cor(result.rgcca$Y[[1]],result.rgcca$Y[[2]])
}
