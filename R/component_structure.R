#A matrix representing common/distinctive processes.
#
#\code{component_structure} generates a matrix representing the common/distinctive components of component loading matrix.
#
#This function generates a matrix indicating which cells in the component loading matrix should be fixed at zeros.
#Thus, it is used together with the function \code{VarSelectComDist}.
#
#@param Jk A vector. Each element of this vector is the number of columns of a data block.
#@param R The number of components.
#@param target A matrix containing 0's and 1's. Its number of columns equals to R, and its number of rows equals to the number of blocks to be integrated. Thus, if the element in
#the first row and first column is 1, then it means that the component belonging to the first block and the first component is selected; if it is 0, then the component is fixed at zeros.
#@return A matrix indicating which cells in the component loading matrix should be fixed at zero.
#@examples
#'target <- matrix(c(1,1,1,0,1,0,0,1,0,1), 2, 5)
#'Jk <- c(144, 44)
#'R <- 5
#'component_structure(Jk, R, target)
#@references Gu, Z., & Van Deun, K. (2016). A variable selection method for simultaneous component based data integration. \emph{Chemometrics and Intelligent Laboratory Systems}, 158, 187-199.

component_structure <- function(Jk, R, target){

  rowTarget <- nrow(target)
  colTarget <- ncol(target)
  if( length(which(target==1 | target==0)) != (rowTarget* colTarget)){
    #something is wrong with the target matrix
    stop("Please enter a proper target matrix with 1 and 0 only. Check help(group_strucure).")
  } else{
    compstructure <- matrix(0, sum(Jk), R)
    for(r in 1:R){
      if (sum(target[,r]) == nrow(target)){
        #common component
        compstructure[, r] <- 1
      } else{
        L <- 1
        for (k in 1:length(Jk)){
          U <- L + Jk[k] - 1
          if(target[k, r] == 1){
            compstructure[L:U, r] <- 1
          } else if(target[k, r] == 0){
            compstructure[L:U, r] <- 0
          }
          L <- U + 1
        }
      }
    }
  }

  return(compstructure)
}
