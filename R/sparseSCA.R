#'Variable selection with Lasso and Group Lasso with a multi-start procedure. 
#'
#'
#'Variable selection with Lasso and Group Lasso penalties to identify component and distinctive components. This algorithm incorporates
#'a multi-start procedure to deal with the possible existence of local minima. 
#'
#'@param DATA A matrix, which contains the concatenated data with the same subjects from multiple blocks.
#'@param Jk A vector containing number of variables in the concatinated data matrix.
#'@param R Number of components (R>=2).
#'@param LASSO A Lasso tuning parameter.
#'@param GROUPLASSO A group Lasso tuning parameter.
#'@param MaxIter The maximum rounds of iterations. It should be a positive integer. The default value is 400.
#'@param NRSTARTS Multi-start procedure: The number of multi-starts. The default value is 20.
#'@param method "datablock" or "component". If \code{method="component"}, the algorithm treats each component across all blocks independently, and thus sparse Group Lasso 
#'is applied per component. If \code{method="datablock"}, the algorithm applies sparse Group Lasso on the entire concatenated data block altogether.
#'If \code{method} is missing, then the \code{"component"} method is used. 
#'
#'@return
#'\item{Pmatrix}{The best estimated component loading matrix (i.e., P), if multi-starts >= 2.}
#'\item{Tmatrix}{The best estimated component score matrix (i.e., T), if multi-starts >= 2.}
#'\item{Lossvec}{A list of vectors containing the loss in each iteration for each multi-start.}
#'
#'@examples
#'\dontrun{
#'DATA1 <- matrix(rnorm(50), nrow=5)
#'DATA2 <- matrix(rnorm(100), nrow=5) 
#'DATA <- cbind(DATA1, DATA2)
#'Jk <- c(10, 20) 
#'R <- 5 
#'LASSO <- 0.2 
#'GROUPLASSO <- 0.4 
#'MaxIter <- 400
#'results <- sparseSCA(DATA, Jk, R, LASSO, GROUPLASSO, 
#'                     MaxIter, NRSTARTS = 10, method = "datablock")
#'
#'results$Pmatrix 
#'}
#'@references
#'Friedman, J., Hastie, T., & Tibshirani, R. (2010). A note on the group lasso and a sparse group lasso. arXiv preprint arXiv:1001.0736.
#'@references
#'Yuan, M., & Lin, Y. (2006). Model selection and estimation in regression with grouped variables. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1), 49-67.
#'@export


sparseSCA <- function(DATA, Jk, R, LASSO, GROUPLASSO, MaxIter, NRSTARTS, method){

  
  if(missing(NRSTARTS)){
    NRSTARTS <- 20
  } 
  
  if(missing(method)){
    method <- "component"
  }
  
  if(R == 1){
    stop("Parameter R = 1 is not allowed.")
  }
  
  Pout3d <- list()
  Tout3d <- list()
  LOSS <- array()
  LOSSvec <- list()
  
  for (n in 1:NRSTARTS){
    
    if(method == "datablock"){
      VarSelectResult <- CDfriedmanV1(DATA, Jk, R, LASSO, GROUPLASSO, MaxIter)
    } else if (method == "component"){
      VarSelectResult <- CDfriedmanV2(DATA, Jk, R, LASSO, GROUPLASSO, MaxIter)
    }
    
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

  return_varselect <- list()
  return_varselect$Pmatrix <- Pout3d[[k]]
  return_varselect$Tmatrix <- Tout3d[[k]]
  return_varselect$Lossvec <- LOSSvec
  attr(return_varselect, "class") <- "sparseSCA"
  return(return_varselect)
}
