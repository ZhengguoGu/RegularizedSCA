#'A K-fold cross-validation procedure when common/distictive processes are
#'known.
#'
#'\code{cv_structuredSCA} helps to find a range of lasso tuning parameters for the
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
#'@param Position A scaler or a vector indicating the column(s) of which the variables will be selected by LASSO in the component loading matrix (P). 
#'@param component_structure A matrix specifing which elements in the component
#'  matrix should be fixed at zeros. see \code{component_structure}.
#'@param MaxIter Maximum number of iterations for this algorithm. The default
#'  value is 400.
#'@param NRSTARTS The number of multistarts for this algorithm. The default
#'  value is 1.
#'@param LassoSequence The range of lasso tuning parameters. The default value
#'  is a sequence of 50 numbers from 0.00000001 to the smallest Lasso tuning parameter
#'  that can make the entire common component(s) to be zeros. Note that by default the 50 numbers are equally spaced on the log scale.
#'@param nfolds Number of folds. If missing, then 10 fold cross-validation will
#'  be performed.
#'@return
#'\item{PRESS}{A vector of predicted residual sum of squares (PRESS) for the sequence of Lasso tuning parameters.}
#'\item{LassoSeqence}{The sequence of Lasso tuning parameters used in cross-validation.}
#'\item{plot}{A plot of mean square errors +/- 1 standard error against Lasso tuning parameters. The plot is plotted agains a log scale of lumbda if \code{LassoSequence} is not defined by users. }
#'\item{LassoRegion}{A region where the suitable lambda can be found, according to the "1 SE rule". }
#'@examples
#'\dontrun{
#'DATA1 <- matrix(rnorm(50), nrow=5)
#'DATA2 <- matrix(rnorm(100), nrow=5) #thus, we assume that DATA1 and DATA2 are with respect to the same 5 subjects here.
#'DATA <- cbind(DATA1, DATA2)
#'Jk <- c(10, 20) #DATA1 has 10 columns, DATA2 20.
#'R <- 4 # assume we want to have 4 components in P matrix.
#'Position <- 1 # assume that we let the variables in the first column to be selected by LASSO in concatenated P matrix.
#'com_str <- component_structure(Jk, R, target) # we can either use the function component_structure() or to specify by ourselves
#'                                                     # here we use the component_structure() function.
#'cv_structuredSCA(DATA, Jk, R, Position, component_structure=com_str, MaxIter = 100, NRSTARTS = 40, LassoSequence = seq(from= 0.002, to=0.1, length.out = 10))
#'# note that since we do now specify nfolds in cv_CDpreKf(), nfolds is set to be 10.
#'}
#'@references Witten, D.M., Tibshirani, R., & Hastie, T. (2009), A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. \emph{Biostatistics}, \emph{10}(3), 515-534.
#'@references Gu, Z., & Van Deun, K. (2016). A variable selection method for simultaneous component based data integration. \emph{Chemometrics and Intelligent Laboratory Systems}, \emph{158}, 187-199.
#'@export
cv_structuredSCA <- function(DATA, Jk, R, Position, component_structure, MaxIter, NRSTARTS, LassoSequence, nfolds){

  DATA <- data.matrix(DATA)
  #this cross-validation function makes use of the CDpre.R.

  if(missing(MaxIter)){
    MaxIter <- 400
  }

  plotlog <- 0 # this is to tell whether the plot is against the log scale of lasso
  if(missing(LassoSequence)){

    results <- CDpre(DATA, Jk, R, Position, component_structure, 0, MaxIter)
    Lassomax <- max(abs(results$Pmatrix[, Position]))
    LassoSequence <- exp(seq(from = log(0.00000001), to = log(Lassomax), length.out = 50))

    plotlog <- 1
  }
  
  if(min(LassoSequence) == 0){
    LassoSequence[which(LassoSequence==min(LassoSequence))] <- 0.00000001
  }
  
  if (min(LassoSequence) < 0) {
    stop("Lasso tuning parameter must be non-negative!")
  }
  if(missing(NRSTARTS)){
    NRSTARTS <- 1
  }

  if(missing(nfolds)){
    nfolds <- 10
  }
  if (nfolds < 2){
    stop("Must be at least 2 folds!")
  }
  PRESS <- array()
  sd_MSE <- array()  #note that this is standard error (not standard deviation, although it's called sd)
  percentRemove <- 1/nfolds
  randmat <- matrix(runif(nrow(DATA) * ncol(DATA)), ncol = ncol(DATA))

  for (l in 1:length(LassoSequence)){
    
    cat(sprintf("\nThe cross-validation procedure might take a while to finish. Please be patient."))
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
        VarSelectResult <- CDpre(DATArm, Jk, R, Position, component_structure, LassoSequence[l], MaxIter)
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
    sd_MSE[l] <- sd(error_x)/sqrt(nfolds)

  }


  upper <- PRESS + sd_MSE
  lower <- PRESS - sd_MSE
  
  if(plotlog == 1){
    df <- data.frame( LassoI = log(LassoSequence), Press = PRESS, Upper = upper, Lower = lower)
    lowestPress <- min(PRESS)
    lowestplus1SE <- lowestPress + sd_MSE[which(min(PRESS) == lowestPress)] #plot 1SE rule region the idea is to fine the region of lasso where according to the 1SE rule the lasso should be in that region. 
    lasso1 <- df$LassoI[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))]
    if(PRESS[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))] - lowestplus1SE > 0 ){
      lasso2 <- df$LassoI[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE))) - 1]
      lregion <- c(lasso2, lasso1)
    } else if (PRESS[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))] - lowestplus1SE < 0 ){
      lasso2 <- df$LassoI[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE))) + 1]
      lregion <- c(lasso1, lasso2)
    } else {
      lasso2 <- lasso1
      lregion <- lasso1
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(x=LassoI,y=Press)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper), width=.1) +
      ggplot2::geom_point(ggplot2::aes(x=LassoI,y=Press)) +
      ggplot2::geom_hline(yintercept = upper[which(PRESS == min(PRESS))], linetype = 3) +
      ggplot2::geom_vline(xintercept = lasso1, 
        linetype = "longdash", col = "red") +
      ggplot2::geom_vline(xintercept = lasso2, 
        linetype = "longdash", col = "red") 
    p <- p + ggplot2::labs(x = "Lasso (log scale)", y="Predicted Mean Squared Errors +/- 1SE")
  } else{
    
    df <- data.frame( LassoI = LassoSequence, Press = PRESS, Upper = upper, Lower = lower)
    
    lowestPress <- min(PRESS)
    lowestplus1SE <- lowestPress + sd_MSE[which(min(PRESS) == lowestPress)] #plot 1SE rule region the idea is to fine the region of lasso where according to the 1SE rule the lasso should be in that region. 
    lasso1 <- df$LassoI[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))]
    if(PRESS[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))] - lowestplus1SE > 0 ){
      lasso2 <- df$LassoI[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE))) - 1]
      lregion <- c(lasso2, lasso1)
    } else if (PRESS[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))] - lowestplus1SE < 0 ){
        if(df$LassoI[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))] == length(LassoSequence)){
          lasso2 <- lasso1
        } else{
          lasso2 <- df$LassoI[which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE))) + 1]
        }
      lregion <- c(lasso1, lasso2)
    } else {
      lasso2 <- lasso1
      lregion <- lasso1
    }
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x=LassoI,y=Press)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower,ymax=Upper), width=.1) +
      ggplot2::geom_point(ggplot2::aes(x=LassoI,y=Press)) +
      ggplot2::geom_hline(yintercept = upper[which(PRESS == min(PRESS))], linetype = 3) +
      ggplot2::geom_vline(xintercept = lasso1, 
        linetype = "longdash", col = "red") +
      ggplot2::geom_vline(xintercept = lasso2, 
        linetype = "longdash", col = "red") 
    p <- p + ggplot2::labs(x = "Lasso", y="Predicted Mean Squared Errors +/- 1SE")
  }
  


  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$LassoSeqence <- LassoSequence
  return_crossvali$plot <- p
  return_crossvali$LassoRegion <- lregion
  return(return_crossvali)

}
