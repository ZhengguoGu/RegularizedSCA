#'

component_strucure <- function(Jk, R, target){

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
