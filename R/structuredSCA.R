#' Variable selection algorithm with a predefined component loading structure
#'
#' Variable selection algorithm when the common/distinctive structure is known a priori.
#' The common component can also be sparse, which is to be estimated by Lasso.
#' The distinctive components are not sparse in the sense that the entire variables in a component (belonging to a certain block) are either all zeros or non-zeros.
#'
#'@param DATA A matrix, which contains the concatenated data with the same subjects from multiple blocks.
#'@param Jk A vector containing number of variables in the concatinated data matrix.
#'@param R Number of components.
#'@param Target A matrix containing 0's and 1's. Its number of columns equals to R, and its number of rows equals to the number of blocks to be integrated. Thus, if the element in
#the first row and first column is 1, then it means that the component belonging to the first block and the first component is selected; if it is 0, then the component is fixed at zeros.
#'@param Position Indicate on which component(s) the Lasso Penalty is imposed. If unspecified, the algorithm assume that the 
#'Lasso penalty is imposed on the common component(s) only. If there is no common component, then Lasso penalty is applied to all components.
#'@param LASSO A Lasso tuning parameter.  
#'@param MaxIter The maximum rounds of iterations. It should be a positive integer. The default value is 400.
#'@param NRSTARTS Multi-start procedure: The number of multi-starts. The default value is 20.
#'
#'@return
#'\item{Pmatrix}{The best estimated component loading matrix (i.e., P), if multi-starts >= 2.}
#'\item{Tmatrix}{The best estimated component score matrix (i.e., T), if multi-starts >= 2.}
#'\item{Lossvec}{A list of vectors containing the loss in each iteration for each multi-start.}
#'
#'@examples
#'\dontrun{
#'DATA1 <- matrix(rnorm(50), nrow=5)
#'DATA2 <- matrix(rnorm(100), nrow=5) #thus, we assume that DATA1 and DATA2 are with respect to the same 5 subjects here.
#'DATA <- cbind(DATA1, DATA2)
#'Jk <- c(10, 20) #DATA1 has 10 columns, DATA2 20.
#'R <- 5 # assume that for some reason, we want to have 5 components in the concatenated P matrix.
#'Target <- matrix(c(1,1,1,0,1,0,0,1,0,1), 2, 5) 
#'LASSO <- 0.2 # assume that we know LASSO=0.2 is a suitable value. 
#'MaxIter <- 400
#'NRSTARTS <- 5
#'structuredSCA(DATA, Jk, R, Target, LASSO = LASSO)
#'}
#'@references
#'Gu, Z., & Van Deun, K. (2016). A variable selection method for simultaneous component based data integration. \emph{Chemometrics and Intelligent Laboratory Systems}, 158, 187-199.
#'@export

structuredSCA <- function(DATA, Jk, R, Target, Position, LASSO, MaxIter, NRSTARTS){
  
  if(missing(NRSTARTS)){
    NRSTARTS <- 20
  } 
  
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  
  GroupStructure <- component_structure(Jk, R, Target)
  
  if(missing(Position)){
    Position <- which(colSums(Target) == nrow(Target))
    
    if(length(Position)==0){
      # no common component
      Position <- 1:R
    }
  }
  
  for (n in 1:NRSTARTS){
    
    VarSelectResult <- CDpre(DATA, Jk, R, Position, GroupStructure, LASSO, MaxIter)
    Pout3d[[n]] <- VarSelectResult$Pmatrix
    Tout3d[[n]] <- VarSelectResult$Tmatrix
    LOSS[n] <- VarSelectResult$Loss
    LOSSvec[[n]] <- VarSelectResult$Lossvec
  }
  
  k <- which(LOSS == min(LOSS))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  PoutBest <- Pout3d[[k]]
  ToutBest <- Tout3d[[k]]
  
  return_varselect <- list()
  return_varselect$Pmatrix <- PoutBest
  return_varselect$Tmatrix <- ToutBest
  return_varselect$Lossvec <- LOSSvec
  return(return_varselect)
}