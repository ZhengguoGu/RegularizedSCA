#' Proportion of variance accounted for (VAF) for each block and each principal component.
#'
#' Proportion of variance accounted for (VAF) is calculated for each block and each column.
#'
#' @param DATA A matrix, which contains the concatenated data with the same subjects from multiple blocks.
#' Note that each row represents a subject.
#' @param Jk A vector containing number of variables in the concatinated data matrix.
#' @param R Number of components (R>=2).
#' @return
#' \item{block}{Proportion of VAF for each block.}
#' \item{component}{Proportion of VAF for each component of each block.}
#'
#' @examples
#' \dontrun{
#' DATA1 <- matrix(rnorm(50), nrow=5)
#' DATA2 <- matrix(rnorm(100), nrow=5) #thus, we assume that DATA1 and DATA2 are with respect to the same 5 subjects here.
#' DATA <- cbind(DATA1, DATA2)
#' Jk <- c(10, 20) #DATA1 has 10 columns, DATA2 20.
#' R <- 5
#' VAF(DATA, Jk, R)
#' }
#' @references
#' Schouteden, M., Van Deun, K., Wilderjans, T. F., & Van Mechelen, I. (2014). Performing DISCO-SCA to search for distinctive and common information in linked data. Behavior research methods, 46(2), 576-587.
#' @references
#' Schouteden, M., Van Deun, K., Pattyn, S., & Van Mechelen, I. (2013). SCA with rotation to distinguish common and distinctive information in linked data. Behavior research methods, 45(3), 822-833.
#'@export
VAF <- function(DATA, Jk, R){

  if(R == 1){
    stop("Parameter R = 1 is not recommended. 
          Please enter an integer higher than 1.")
  }
  
  DATA <- data.matrix(DATA)
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
