#'A K-fold cross-validation procedure when common/distictive processes are unknown.
#'
#'\code{cv_CDfriedman} helps to find a range of Lasso and Group Lasso tuning parameters for the common component so as to generate sparse common component.
#'
#'This function search through a range of Lasso and Group Lasso tuning parameters for identifying common and distinctive components
#'
#'@param DATA The concatenated data block, with rows represending subjects.
#'@param Jk A vector. Each element of this vector is the number of columns of a data block.
#'@param R The number of components.
#'@param MaxIter Maximum number of iterations for this algorithm. The default value is 400.
#'@param NRSTARTS The number of multistarts for this algorithm. The default value is 5.
#'@param LassoSequence The range of Lasso tuning parameters. The default value is a sequence of 5 numbers from 0.001
#'to the smallest Lasso tuning parameter that can make all the components to be zeros.
#'@param GLassoSequence The range of Group Lasso tuning parameters. The default value is a sequence of 5 numbers from 0.001
#'to the smallest Group Lasso tuning parameter that can make all the components to be zeros.
#'@param nfolds Number of folds. If missing, then 10 fold cross-validation will be performed.
#'@return
#'\item{PRESS}{A matrix of predicted residual sum of squares (PRESS) for the sequences of Lasso and Group Lasso tuning parameters.}
#'\item{LassoSequence}{The sequence of Lasso tuning parameters used in cross-validation.}
#'\item{GLassoSequence}{The sequence of Group Lasso tuning parameters used in cross-validation.}
#'\item{plot}{A plot of PRESS against Lasso and Group Lasso tuning parameters. }
#'\item{plotSE}{A plot of PRESS +/- 1 standard error against Lasso and Group Lasso tuning parameters. }
#'
#'@examples
#'cv.CDpre(DATA, Jk, R, CommPosition, component_structure, MaxIter = 100, NRSTARTS = 40, LassoSequence = seq(from= 0.002, to=0.1, length.out = 10))
#'@references Witten, D.M., Tibshirani, R., & Hastie, T. (2009), A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. \emph{Biostatistics}, \emph{10}(3), 515-534.
#'@references
#'Friedman, J., Hastie, T., & Tibshirani, R. (2010). A note on the group lasso and a sparse group lasso. arXiv preprint arXiv:1001.0736.
#'@references
#'Yuan, M., & Lin, Y. (2006). Model selection and estimation in regression with grouped variables. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1), 49-67.
cv_CDfriedman <- function(DATA, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds){

  if(missing(LassoSequence) | missing(GLassoSequence)){

      results <- maxLGlasso(DATA, Jk, R)
      GLassomax <- results$Glasso
      Lassomax <- results$Lasso

    if(missing(LassoSequence)){
      LassoSequence <- seq(from = 0.001, to = Lassomax, length.out = 5)
    }

    if(missing(GLassoSequence)){
      GLassoSequence <- seq(from = 0.001, to = GLassomax, length.out = 5)
      }
  }

  if(missing(MaxIter)){
    MaxIter <- 400
  }

  if(missing(NRSTARTS)){
    NRSTARTS <- 5
  }

  if(missing(nfolds)){
    nfolds <- 10
  }
  if (nfolds < 2)
    stop("Must be at least 2 folds")

  PRESS <- matrix(0, length(LassoSequence), length(GLassoSequence))
  se_MSE <- matrix(0, length(LassoSequence), length(GLassoSequence))
  percentRemove <- 1/nfolds
  randmat <- matrix(runif(nrow(DATA) * ncol(DATA)), ncol = ncol(DATA))

  for (g in 1: length(GLassoSequence)){

    for (l in 1:length(LassoSequence)){

      cat(sprintf("\nThe cross-validation procedure might take while to finish. Please be patient."))
      cat(sprintf("\nGroup Lasso: %s", GLassoSequence[g]))
      cat(sprintf("\nLasso: %s", LassoSequence[l]))

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
          VarSelectResult <- CDfriedman(DATArm, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
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


      PRESS[l,g] <- sum(error_x)/nfolds
      se_MSE[l,g] <- sd(error_x)/sqrt(nfolds)
    }
  }

  vec_PRESS <- c(PRESS)
  vec_se <- c(se_MSE)
  #lasso_label <- rep(LassoSequence, length(GLassoSequence))
  lasso_index <- rep(1:length(LassoSequence), length(GLassoSequence))
  #Glasso_label <- rep(GLassoSequence, each=length(LassoSequence))
  Glasso_index <- rep(1:length(GLassoSequence), each=length(LassoSequence))


  plot(1:length(lasso_index), vec_PRESS, xlab = "" , ylab = "", axes=FALSE)
  Lseq <- 1:length(lasso_index)
  abline(v=Lseq[which(lasso_index==max(lasso_index))], lty=2)
  axis(2, round(seq(from=min(vec_PRESS), to=max(vec_PRESS), length.out = 10), digits = 3), las=1)
  axis(1,at=1:length(lasso_index), labels = lasso_index, line=0)
  mtext("Lasso",1,line=0,at=-1)
  axis(1,at=1:length(Glasso_index), labels = Glasso_index, line=2)
  mtext("G-Lasso",1,line=2,at=-1)
  pic1 <- recordPlot()

  y_min <- min(vec_PRESS-vec_se)
  y_max <- max(vec_PRESS+vec_se)
  plot(1:length(lasso_index), vec_PRESS, xlab = "" , ylab = "", axes=FALSE, xaxt='n', yaxt="n")
  abline(v=Lseq[which(lasso_index==max(lasso_index))], lty=2)
  arrows(1:length(lasso_index), vec_PRESS-vec_se, 1:length(lasso_index), vec_PRESS+vec_se, length=0.05, angle=90, code=3)
  axis(2, round(seq(from=y_min, to=y_max, length.out = 10), digits = 3), las=1)
  axis(1,at=1:length(lasso_index), labels = lasso_index, line=0)
  mtext("Lasso",1,line=0,at=-1)
  axis(1,at=1:length(Glasso_index), labels = Glasso_index, line=2)
  mtext("G-Lasso",1,line=2,at=-1)
  pic2 <- recordPlot()

  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$LassoSeqence <- LassoSequence
  return_crossvali$GLassoSeqence <- GLassoSequence
  return_crossvali$plot <- pic1
  return_crossvali$plotSE <- pic2
  return(return_crossvali)

}

