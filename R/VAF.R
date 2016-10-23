# proportion of variance accounted for

VAF <- function(DATA, Jk, R){

  SVD_DATA <- svd(DATA, R, R)
  Tmat <- SVD_DATA$u
  Pmat <- SVD_DATA$v %*% diag(SVD_DATA$d[1:R])

  VAF_results_block <- array(NA, length(Jk))
  VAF_results_component <- matrix(NA, length(Jk), R)

  L <- 1
  for (i in 1:length(Jk)){

    U <- L + Jk[i] - 1
    X_hat <- Tmat%*%t(Pmat[L:U, ])
    DATA_k <- DATA[, L:U]

    VAF_results_block[i] <- 1 - sum((X_hat - DATA_k)^2)/sum(DATA_k^2)

    Pmat_k <- Pmat[L:U, ]
    for (r in 1: R){
      x_pwise_hat <- Tmat[, r] %*% t(Pmat_k[, r])
      VAF_results_component[i, r] <- sum(x_pwise_hat^2)/sum(DATA_k^2)
    }
    L <- U + 1
  }
  return_VAF <- list()
  return_VAF$block <- VAF_results_block
  return_VAF$component <- VAF_results_component

  return(return_VAF)
}
