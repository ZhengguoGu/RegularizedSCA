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
#'\item{plot}{A plot of PRESS against Lasso and Group Lasso tuning parameters. Note that on the x axis are the index numbers of
#'Lasso and Group Lasso tuning parameters, and to find their corresponding values, please make use of \code{lasso_index} and \code{Glasso_index}}
#'\item{plotSE}{A plot of PRESS +/- 1 standard error against Lasso and Group Lasso tuning parameters. Note that on the x axis are the index numbers of
#'Lasso and Group Lasso tuning parameters, and to find their corresponding values, please make use of \code{lasso_index} and \code{Glasso_index}}
#'\item{Lasso_values}{The sequence of Lasso tuning parameters used for cross-validation. For example, suppose from the plot we found that the index number
#'for Lasso is \code{6}, its corresponding Lasso tuning parameter is \code{Lasso_values[6]}.}
#'\item{Glasso_values}{The sequence of Group Lasso tuning parameters used for cross-validation. For example, suppose from the plot we found that the index number
#'for Group Lasso is \code{6}, its corresponding Group Lasso tuning parameter is \code{Glasso_values[6]}.}
#'
#'@examples
#'\dontrun{
#'DATA1 <- matrix(rnorm(50), nrow=5)
#'DATA2 <- matrix(rnorm(100), nrow=5) #thus, we assume that DATA1 and DATA2 are with respect to the same 5 subjects here.
#'DATA <- cbind(DATA1, DATA2)
#'Jk <- c(10, 20) #DATA1 has 10 columns, DATA2 20.
#'# assume that we do not know which values to choose for the Lasso/Group Lasso tuning parameter, then no need to specify them.
#'cv_CDfriedman(DATA, Jk, R=5, MaxIter = 100, NRSTARTS = 40, nfolds=10)
#'}
#'@references Witten, D.M., Tibshirani, R., & Hastie, T. (2009), A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. \emph{Biostatistics}, \emph{10}(3), 515-534.
#'@references
#'Friedman, J., Hastie, T., & Tibshirani, R. (2010). A note on the group lasso and a sparse group lasso. arXiv preprint arXiv:1001.0736.
#'@references
#'Yuan, M., & Lin, Y. (2006). Model selection and estimation in regression with grouped variables. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1), 49-67.
#'@export
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

  ii <- 0
  while(ii != 1){ #this procedure is to make sure that the training sample do not have an entire row/column of NA's

    randmat <- matrix(runif(nrow(DATA) * ncol(DATA)), ncol = ncol(DATA))
    jj <- 0
    for (i in 1:nfolds){
      ToRemove <- ((i - 1) * percentRemove < randmat) & (randmat < i * percentRemove) # this idea is from PMA package
      DATArm <- DATA
      DATArm[ToRemove] <- NA

      if(ncol(DATArm) %in% rowSums(is.na(DATArm))){
        break
      } else if(nrow(DATArm) %in% colSums(is.na(DATArm))){
        break
      }else{
        jj <- jj+1
      }
    }

    if(jj == nfolds){
      ii <- 1
    }
  }

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
  axis(2, round(seq(from=min(vec_PRESS), to=max(vec_PRESS), length.out = 10), digits = 2), las=1)
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
  axis(2, round(seq(from=y_min, to=y_max, length.out = 10), digits = 2), las=1)
  axis(1,at=1:length(lasso_index), labels = lasso_index, line=0)
  mtext("Lasso",1,line=0,at=-1)
  axis(1,at=1:length(Glasso_index), labels = Glasso_index, line=2)
  mtext("G-Lasso",1,line=2,at=-1)
  pic2 <- recordPlot()

  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$plot <- pic1
  return_crossvali$plotSE <- pic2
  return_crossvali$Lasso_values <- LassoSequence
  return_crossvali$Glasso_values <- GLassoSequence
  return(return_crossvali)

}

