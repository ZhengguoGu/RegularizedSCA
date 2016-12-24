crossvali <- function(DATA, Jk, R, CommPosition, GroupStructure, MaxIter, NRSTARTS, LassoSequence){

  #this cross-validation function makes use of the VarSelectComDistPre.R.
  if(missing(LassoSequence)){

    results <- VarSelectComDistPre(DATA, Jk, R, CommPosition, GroupStructure, 0, MaxIter)
    Lassomax <- max(abs(results$Pmatrix[, CommPosition]))
    LassoSequence <- seq(from = 0, to = Lassomax, length.out = 10)

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
          VarSelectResult <- VarSelectComDistPre(DATA_x, Jk, R, CommPosition, GroupStructure, LassoSequence[l], MaxIter)
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

  plot(LassoSequence[-1], PRESS[-1], xlab = 'Lasso tuning parameter', ylab = 'Predicted Residual Sum of Squares')
  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$LassoSeqence <- LassoSequence
  return(return_crossvali)
}
