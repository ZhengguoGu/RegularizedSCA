#' Tucker congruence
#'
#' \code{TuckerCoef} calculate Tucker's coefficient of congruence between columns but after accounting for permutational freedom and reflections
#'
#' @param MatrixA A matrix
#' @param MatrixB A matrix, which is to be compared to MatrixA
#' @return
#' \item{perm}{the permutation order.}
#' \item{tucker_value}{the Tucker coefficient.}
#' \item{tucker_vector}{the Tucker vector.}
#' @examples
#' maxtrix1 <- matrix(rnorm(50), nrow=5)
#' maxtrix2 <- matrix(rnorm(50), nrow=5)
#' TuckerCoef(maxtrix1, maxtrix2)
#' @references 
#' Lorenzo-Seva, U., & Ten Berge, J. M. (2006). Tucker's congruence coefficient as a meaningful index of factor similarity. \emph{Methodology}, \emph{2}(2), 57-64.
#'@export
TuckerCoef <- function(MatrixA, MatrixB){
  
  nrow_data <- dim(MatrixA)[1]
  ncol_data <- dim(MatrixA)[2]
  INDIC_Mat <- gtools::permutations(ncol_data, ncol_data)
  ncol_INDIC <- dim(INDIC_Mat)[1]
  TUCK <- array(NA, dim = c(ncol_INDIC, ncol_data))
  tucker_values <- array()
  tuckerr <- array()
  for(i in 1: ncol_INDIC) {
    MatrixB_perm <- MatrixB[, INDIC_Mat[i,]]
    teller <- 1

    for (r in 1: ncol_data){
      vec1 <- MatrixA[, r]
      vec2 <- MatrixB_perm[, r]
      cp <- t(vec1) %*% vec2
      var1 <- t(vec1) %*% vec1
      var2 <- t(vec2) %*% vec2

      if (var1 > 0 & var2 > 0){
        tuckerr[teller] <- psych::tr(cp)/sqrt(psych::tr(var1)*psych::tr(var2))
        teller <- teller + 1
      } else if (var2 == 0){
        tuckerr[teller] <- 0
        teller <- teller + 1
      }
    }

    tucker_values[i] <- mean(abs(tuckerr))
    TUCK[i,] <- tuckerr
  }

  k <- which(tucker_values == max(tucker_values))
  k <- k[1]

  perm <- INDIC_Mat[k,]
  tucker_value <- max(tucker_values)
  tucker_vector <- TUCK[k, ]

  return_tucker <- list()
  return_tucker$perm <- perm
  return_tucker$tucker_value <- tucker_value
  return_tucker$tucker_vector <- tucker_vector
  return(return_tucker)
}
