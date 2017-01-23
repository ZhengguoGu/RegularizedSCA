#'A K-fold cross-validation procedure when common/distictive processes are
#'known.
#'
#'\code{cv_CDpreKf} helps to find a range of lasso tuning parameters for the
#'common component so as to generate sparse common component.
#'
#'This function search through a range of lasso tuning parameters for the common
#'component, while keeping distinctive components fixed (- that is, the zeros in
#'the distinctive components are fixed). This function may be of help if user
#'wants to obtain some sparsness in the common component.
#'
#'@param DATA The concatenated data block, with rows represending subjects.
#'@param Jk A vector. Each element of this vector is the number of columns of a
#'  data block.
#'@param R The number of components.
#'@param CommPosition A number (vector) indicating which component(s) is the
#'  common component(s).
#'@param component_structure A matrix specifing which elements in the component
#'  matrix should be fixed at zeros. see \code{component_structure}.
#'@param MaxIter Maximum number of iterations for this algorithm. The default
#'  value is 400.
#'@param NRSTARTS The number of multistarts for this algorithm. The default
#'  value is 10.
#'@param LassoSequence The range of lasso tuning parameters. The default value
#'  is a sequence of 10 numbers from 0 to the smallest Lasso tuning parameter
#'  that can make the entire common component(s) to be zeros.
#'@param nfolds Number of folds. If missing, then 10 fold cross-validation will
#'  be performed.

#'@return
#'\item{PRESS}{A vector of predicted residual sum of squares (PRESS) for the sequence of Lasso tuning parameters.}
#'\item{LassoSeqence}{The sequence of Lasso tuning parameters used in cross-validation.}
#'\item{plot}{A plot of mean square errors against Lasso tuning parameters.}
#'\item{plotSD}{A plot of mean square errors +/- 1SD against Lasso tuning parameters.}
#'@examples
#'cv.CDpre(DATA, Jk, R, CommPosition, component_structure, MaxIter = 100, NRSTARTS = 40, LassoSequence = seq(from= 0.002, to=0.1, length.out = 10))
#'@references Witten, D.M., Tibshirani, R., & Hastie, T. (2009), A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. \emph{Biostatistics}, \emph{10}(3), 515-534.
#'@references Gu, Z., & Van Deun, K. (2016). A variable selection method for simultaneous component based data integration. \emph{Chemometrics and Intelligent Laboratory Systems}, \emph{158}, 187-199.
cv_CDpreKf <- function(DATA, Jk, R, CommPosition, component_structure, MaxIter, NRSTARTS, LassoSequence, nfolds){

  #this cross-validation function makes use of the CDpre.R.

  if(missing(MaxIter)){
    MaxIter <- 400
  }

  if(missing(LassoSequence)){

    results <- CDpre(DATA, Jk, R, CommPosition, component_structure, 0, MaxIter)
    Lassomax <- max(abs(results$Pmatrix[, CommPosition]))
    LassoSequence <- seq(from = 0, to = Lassomax, length.out = 10)

  }

  if(missing(NRSTARTS)){
    NRSTARTS <- 10
  }

  if(missing(nfolds)){
    nfolds <- 10
  }
  if (nfolds < 2)
    stop("Must be at least 2 folds")

  PRESS <- array()
  sd_MSE <- array()
  percentRemove <- 1/nfolds
  randmat <- matrix(runif(nrow(DATA) * ncol(DATA)), ncol = ncol(DATA))

  for (l in 1:length(LassoSequence)){
    error_x <- array()

    for (i in 1:nfolds){
      ToRemove <- ((i - 1) * percentRemove < randmat) & (randmat < i * percentRemove) # this idea is from PMA package
      DATArm <- DATA
      DATArm[ToRemove] <- NA

      for(c in 1:ncol(DATA)){
        indexc <- !is.na(DATArm[, c])
        DATArm[, c][!indexc] <- mean(DATArm[, c][indexc]) #missing values are replaced by column means
      }

      Pout3d <- list()
      Tout3d <- list()
      LOSS <- array()
      for (n in 1:NRSTARTS){
        VarSelectResult <- CDpre(DATArm, Jk, R, CommPosition, component_structure, LassoSequence[l], MaxIter)
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
      error_x[i] <- sum((DATA[ToRemove] - DATA_hat[ToRemove])^2)

    }

    PRESS[l] <- sum(error_x)/nfolds
    sd_MSE[l] <- sd(error_x)

  }



  plot(LassoSequence, PRESS, xlab = 'Lasso tuning parameter', ylab = 'Mean Square Error', type = "b")
  pic1 <- recordPlot()

  y_min <- min(PRESS-sd_MSE)
  y_max <- max(PRESS+sd_MSE)
  plot(LassoSequence, PRESS, xlab = 'Lasso tuning parameter', ylab = 'Mean Square Error +/- 1SD', ylim = c(y_min, y_max), type = "b")
  arrows(LassoSequence, PRESS-sd_MSE, LassoSequence, PRESS+sd_MSE, length=0.05, angle=90, code=3)
  pic2 <- recordPlot()

  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$LassoSeqence <- LassoSequence
  return_crossvali$plot <- pic1
  return_crossvali$plotSD <- pic2
  return(return_crossvali)

}
