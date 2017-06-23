###################################################
#  a function to construct matrix of reflections
#  author: Katrijn Van Duen
#  implemented by Zhengguo Gu
###################################################

reflexmat <- function(m){

  mat <- rep(1, m)

  for(i in 1:(m-1)){

    B <- utils::combn(1:m, i)

      for(j in 1:dim(B)[2]){

        v <- rep(1, m)
        v[t(B[, j])] <- -1

        mat <- rbind(mat, v)
      }
  }

  return(mat)
}
