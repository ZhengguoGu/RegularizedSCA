#' Undo shrinkage. 
#' 
#' \code{undoShrinkage} re-estimates the component loading matrix (P) while keeping the 0 loadings fixed 
#' so as to remove the shrinkage due to Lasso and Group Lasso. 
#' 
#' @param DATA The concatenated data block, with rows representing subjects
#' @param R The number of components.
#' @param Phat The estimated component loading matrix by means of, for example, \code{sparseSCA()}. 
#' @param MAXITER The maximum rounds of iterations. It should be a positive integer. The default value is 400.
#' 
#' @return 
#' \item{Pmatrix}{The re-estimated component loading matrix after the shrinkage has been removed.}
#' \item{Tmatrix}{The corresponding estimated component score matrix.}
#' \item{Lossvec}{A vector of loss.}
#' 
#'@references
#'Gu, Z., & Van Deun, K. (2016). A variable selection method for simultaneous component based data integration. \emph{Chemometrics and Intelligent Laboratory Systems}, 158, 187-199.
#'@export


undoShrinkage <- function(DATA, R, Phat, MAXITER){
  
  DATA <- data.matrix(DATA)
  I_Data <- dim(DATA)[1]
  sumJk <- dim(DATA)[2]
  eps <- 10^(-12)
  
  if(missing(MAXITER)){
    MAXITER <- 400
  }
  
  #initialize P, Lossc
  P <- matrix(stats::rnorm(sumJk * R), nrow = sumJk, ncol = R)
  P[Phat == 0]<-0
  Pt <- t(P)
  
  residual <- sum(DATA^2)
  Lossc <- residual
  Lossvec <- array()
  
  conv <- 0
  iter <- 1
  while (conv == 0){
    
    #update T
    A <- Pt %*% t(DATA)
    SVD_DATA <- svd(A, R, R)
    Tmat <- SVD_DATA$v %*% t(SVD_DATA$u)
    
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu <- residual
    
    #update P
    
    P <- t(DATA) %*% Tmat
    P[Phat == 0] <- 0  # this is to keep the zero structures
    Pt <- t(P)
    
    residual <- sum((DATA - Tmat %*% Pt)^2)
    Lossu2 <- residual
    
    #check convergence
    if (abs(Lossc-Lossu)< 10^(-9)) {
      Loss <- Lossu
      residual <- residual
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    else if (iter > MAXITER){
      Loss <- Lossu
      residual <- residual
      P[abs(P) <= 2 * eps] <- 0
      conv <- 1
    }
    
    Lossvec[iter] <- Lossu
    iter <- iter + 1
    Lossc <- Lossu2
    
  }
  
  return_undoShrink <- list()
  return_undoShrink$Pmatrix <- P
  return_undoShrink$Tmatrix <- Tmat
  #return_undoShrink$Loss <- Loss
  return_undoShrink$Lossvec <- Lossvec
  attr(return_undoShrink, "class") <- "undoS"
  
  return(return_undoShrink)
}
