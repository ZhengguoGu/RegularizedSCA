#'A cross-validation procedure when common/distictive processes are known.
#'
#'\code{crossvali} helps to find a range of lasso tuning parameters for the common component so as to generate sparse common component.
#'
#'This function search through a range of lasso tuning parameters for the common component, while keeping
#'distinctive components fixed (- that is, the zeros in the distinctive components are fixed). This function
#'may be of help if user wants to obtain some sparsness in the common component.
#'
#'@param DATA The concatenated data block, with rows represending subjects.
#'@param Jk A vector. Each element of this vector is the number of columns of a data block.
#'@param R The number of components.
#'@param CommPosition A number (vector) indicating which component(s) is the common component(s).
#'@param component_structure A matrix specifing which elements in the component matrix should be fixed at zeros.
#'see \code{component_structure}.
#'@param MaxIter Maximum number of iterations for this algorithm. The default value is 400.
#'@param NRSTARTS The number of multistarts for this algorithm. The default value is 20.
#'@param LassoSequence The range of lasso tuning parameters. The default value is a sequence of 10 numbers from 0
#'to the smallest Lasso tuning parameter that can make the entire common component(s) to be zeros.
#'@return
#'\item{PRESS}{A vector of predicted residual sum of squares (PRESS) for the sequence of Lasso tuning parameters.}
#'\item{LassoSeqence}{The sequence of Lasso tuning parameters used in cross-validation.}
#'\item{plot}{A plot of PRESS against Lasso tuning parameters. }
#'
#'@examples
#'crossvali(DATA, Jk, R, CommPosition, component_structure, MaxIter = 100, NRSTARTS = 40, LassoSequence = seq(from= 0.002, to=0.1, length.out = 10))
#'@references Bro, R., Kjeldahl, K., Smilde, A. K., & Kiers, H. A. L. (2008). Cross-validation of component models: a critical look at current methods. Analytical and bioanalytical chemistry, 390(5), 1241-1251.
#'@references Gu, Z., & Van Deun, K. (2016). A variable selection method for simultaneous component based data integration. \emph{Chemometrics and Intelligent Laboratory Systems}, 158, 187-199.
crossvali <- function(DATA, Jk, R, CommPosition, component_structure, MaxIter, NRSTARTS, LassoSequence){

  #this cross-validation function makes use of the CDpre.R.
  if(missing(LassoSequence)){

    results <- CDpre(DATA, Jk, R, CommPosition, component_structure, 0, MaxIter)
    Lassomax <- max(abs(results$Pmatrix[, CommPosition]))
    LassoSequence <- seq(from = 0, to = Lassomax, length.out = 10)

  }

  if(missing(MaxIter)){
    MaxIter <- 400
  }

  if(missing(NRSTARTS)){
    NRSTARTS <- 20
  }

  PRESS <- array()
  for (l in 1:length(LassoSequence)){
    error_x <- 0
    for (i in nrow(DATA)){
      for (j in sum(Jk)){

        Pout3d <- list()
        Tout3d <- list()
        LOSS <- array()
        DATA_x <- DATA
        x <- DATA[i,j]
        DATA_x[i, j] <- 0 #remove the data point

        for (n in 1:NRSTARTS){
          VarSelectResult <- CDpre(DATA_x, Jk, R, CommPosition, component_structure, LassoSequence[l], MaxIter)
          Pout3d[[n]] <- VarSelectResult$Pmatrix
          Tout3d[[n]] <- VarSelectResult$Tmatrix
          LOSS[n] <- VarSelectResult$Loss
        }
        k <- which(LOSS == min(LOSS))
        if (length(k)>1){
          pos <- sample(1:length(k), 1)
          k <- k[pos]
          }
        PoutBest <- Pout3d[[k]]
        ToutBest <- Tout3d[[k]]

        DATA_hat <- ToutBest%*%t(PoutBest)
        error_x <- error_x + (x - DATA_hat[i, j])^2
      }
    }

    PRESS[l] <- error_x/(nrow(DATA)*sum(Jk))
  }

  pic <- plot(LassoSequence, PRESS, xlab = 'Lasso tuning parameter', ylab = 'Predicted Residual Sum of Squares')
  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$LassoSeqence <- LassoSequence
  return_crossvali$plot <- pic
  return(return_crossvali)
}
