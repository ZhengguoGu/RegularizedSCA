#'A K-fold cross-validation procedure when common/distinctive processes are unknown with Lasso and Group Lasso penalties.
#'
#'\code{cv_sparseSCA} helps to find a range of Lasso and Group Lasso tuning parameters for the common component so as to generate sparse common component.
#'
#'This function searches through a range of Lasso and Group Lasso tuning parameters for identifying common and distinctive components
#'
#'@param DATA The concatenated data block, with rows representing subjects.
#'@param Jk A vector. Each element of this vector is the number of columns of a data block.
#'@param R The number of components (R>=2).
#'@param MaxIter Maximum number of iterations for this algorithm. The default value is 400.
#'@param NRSTARTS The number of multistarts for this algorithm. The default value is 1.
#'@param LassoSequence The range of Lasso tuning parameters. The default value is a sequence of 20 numbers from 0.00000001
#'to the smallest Lasso tuning parameter value that makes all the component loadings equal to zero. Note that by default the 20 numbers are equally spaced on the log scale. 
#'Furthermore, if \code{GLassoSequence} contains only one number, then by default \code{LassoSequence} is a sequence of 50 values.
#'@param GLassoSequence The range of Group Lasso tuning parameters. The default value is a sequence of 10 numbers from 0.00000001
#'to the smallest Group Lasso tuning parameter value that makes all the component loadings equal to zero. Note that by default the 10 numbers are equally spaced (but not on the log scale). 
#'Note that if \code{LassoSequence} contains only one number, then by default \code{GLassoSequence} is a sequence of 50 values.
#'@param nfolds Number of folds. If missing, then 10 fold cross-validation will be performed.
#'@param method "datablock" or "component". These are two options with respect to the grouping of the loadings as used in the Group Lasso penalty. 
#'If \code{method="component"}, the block-grouping of the coefficients is applied per component separately. If \code{method = "datablock"}, the grouping
#'is applied on the concatenated data block, with loadings of all components together. If \code{method} is missing, then the "component" method is used 
#'by default.  
#'
#'@return
#'\item{PRESS}{A matrix of predicted residual sum of squares (PRESS) for the sequences of Lasso and Group Lasso tuning parameters.}
#'\item{SE_MSE}{A matrix of standard errors for \code{PRESS}.}
#'\item{Press1SE}{The lowest PRESS + 1SE.}
#'\item{VarSelected}{A matrix of number of variables selected for the sequences of Lasso and Group Lasso tuning parameters.}
#'\item{plot}{A plot of PRESS +/- 1 standard error against Lasso and Group Lasso tuning parameters, with the vertical dotted black line indicating the lowest
#'PRESS+1SE. Note that on the x axis (abscissa) are Lasso tuning parameter values. The Group Lasso tuning parameter values are shown on the top of the graph, and the values shown are index numbers:
#'G1, for example, indicates the first value in the \code{GLassoSequence}.
#'In case both the Lasso sequence and the Group Lasso sequence contain more than 2 elements, there will be an extra plot, which is 
#'the number of non-zero component loadings against Lasso and Group Lasso tuning parameters. In this case \code{plot} is a list of two plots.
#'To find their corresponding values, please make use of \code{Lasso_values} and \code{Glasso_values}. The vertical red dashed lines indicate a proper region for Lasso tuning parameters
#'given a certain Group Lasso tuning parameter. When there is only one vertical red dashed line, a proper region for Lasso tuning parameters is not available: the red dashed line indicates 
#'the Lasso tuning parameter leading to the PRESS that is closest to the smallest PRESS+1SE.} 
#'\item{Lasso_values}{The sequence of Lasso tuning parameters used for cross-validation. Users may also consult \code{Lambdaregion} (explained below).}
#'\item{Glasso_values}{The sequence of Group Lasso tuning parameters used for cross-validation. For example, suppose from the plot we found that the index number
#'for Group Lasso is \code{6}, its corresponding Group Lasso tuning parameter is \code{Glasso_values[6]}.}
#'\item{Lambdaregion}{A region of proper tuning parameter values for Lasso, given a certain value for Group Lasso. This means that, for example, if 5 Group Lasso tuning parameter values have been considered, \code{Lambdaregion} is a 5 by 2 matrix.}
#'\item{RecommendedLambda}{A pair (or sometimes a few pairs) of Lasso and Group Lasso tuning parameters that lead to a model with PRESS closest to the lowest PRESS + 1SE.}
#'\item{plotlog}{An index number for function \code{plot()}, which is not useful for users.}
#'@examples
#'\dontrun{
#'DATA1 <- matrix(rnorm(50), nrow=5)
#'DATA2 <- matrix(rnorm(100), nrow=5)  
#'DATA <- cbind(DATA1, DATA2)
#'Jk <- c(10, 20) 
#'cv_sparseSCA(DATA, Jk, R=5, MaxIter = 100, NRSTARTS = 40, nfolds=10)
#'}
#'@references Witten, D.M., Tibshirani, R., & Hastie, T. (2009), A penalized matrix decomposition, with applications to sparse principal components and canonical correlation analysis. \emph{Biostatistics}, \emph{10}(3), 515-534.
#'@references
#'Friedman, J., Hastie, T., & Tibshirani, R. (2010). A note on the group lasso and a sparse group lasso. arXiv preprint arXiv:1001.0736.
#'@references
#'Yuan, M., & Lin, Y. (2006). Model selection and estimation in regression with grouped variables. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 68(1), 49-67.
#'@export
cv_sparseSCA <- function(DATA, Jk, R, MaxIter, NRSTARTS, LassoSequence, GLassoSequence, nfolds, method){

  DATA <- data.matrix(DATA)
  plotlog <- 0

  if(R == 1){
    stop("Parameter R = 1 is not allowed.")
  }
  
  if(missing(LassoSequence) | missing(GLassoSequence)){

      results <- maxLGlasso(DATA, Jk, R)
      GLassomax <- results$Glasso
      Lassomax <- results$Lasso

      
    if(missing(LassoSequence) & missing(GLassoSequence)){
      LassoSequence <- exp(seq(from = log(0.00000001), to = log(Lassomax), length.out = 20))
      GLassoSequence <- seq(from = 0.00000001, to = GLassomax, length.out = 10)  #note that Glasso is not on the log scale, because it is not helpful.
      plotlog <- 1
    } else if(missing(LassoSequence) & (length(GLassoSequence) == 1)){
        LassoSequence <- exp(seq(from = log(0.00000001), to = log(Lassomax), length.out = 50))
        plotlog <- 1
    } else if(missing(GLassoSequence) & (length(LassoSequence) == 1)){
        GLassoSequence <- seq(from = 0.00000001, to = GLassomax, length.out = 50)
    } else if(missing(GLassoSequence) & (length(LassoSequence)>1)){
        GLassoSequence <- seq(from = 0.00000001, to = GLassomax, length.out = 10)
    } else if(missing(LassoSequence) & (length(GLassoSequence)>1)){
        LassoSequence <- exp(seq(from = log(0.00000001), to = log(Lassomax), length.out = 20))
        plotlog <- 1
    }
    
    
  }

  if(length(Jk) <= 1){
    stop("This is an algorithm for data integration! There should be at least 2 data blocks!")
  }
  
  if (min(GLassoSequence) < 0) {
    stop("Group Lasso tuning parameter must be non-negative!")
  }
  
  if (min(LassoSequence) < 0) {
    stop("Lasso tuning parameter must be non-negative!")
  }
  if(missing(MaxIter)){
    MaxIter <- 400
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
  
  if(missing(method)){
    method <- "component"
  }
  
  
  PRESS <- matrix(0, length(LassoSequence), length(GLassoSequence))
  se_MSE <- matrix(0, length(LassoSequence), length(GLassoSequence))
  varselected <- matrix(0, length(LassoSequence), length(GLassoSequence))
  percentRemove <- 1/nfolds

  ii <- 0
  while(ii != 1){ #this procedure is to make sure that the training sample do not have an entire row/column of NA's

    randmat <- matrix(stats::runif(nrow(DATA) * ncol(DATA)), ncol = ncol(DATA))
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

      cat(sprintf("\nThe cross-validation procedure might take a while to finish. Please be patient."))
      cat(sprintf("\nGroup Lasso: %s", GLassoSequence[g]))
      cat(sprintf("\nLasso: %s", LassoSequence[l]))

      if(method == "datablock"){
        Forvarselected <- CDfriedmanV1(DATA, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
      }else if (method == "component"){
        Forvarselected <- CDfriedmanV2(DATA, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
      }
      varselected[l,g] <- sum(Forvarselected$Pmatrix != 0)  #how many variables in P have been selected?
      
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
          if(method == "datablock"){
            VarSelectResult <- CDfriedmanV1(DATArm, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
          }else if (method == "component"){
            VarSelectResult <- CDfriedmanV2(DATArm, Jk, R, LassoSequence[l], GLassoSequence[g], MaxIter)
          }
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
      se_MSE[l,g] <- stats::sd(error_x)/sqrt(nfolds)
    }
  }

  vec_PRESS <- c(PRESS)
  vec_se <- c(se_MSE)
  vec_varsel <- c(varselected)
  upper <- vec_PRESS + vec_se
  lower <- vec_PRESS - vec_se 
  #lasso_index0 <- rep(1:length(LassoSequence), length(GLassoSequence))
  #Glasso_index0 <- rep(1:length(GLassoSequence), each=length(LassoSequence))
  
  #lasso_index <- paste("L", lasso_index0)
  #Glasso_index<- factor(paste("G", Glasso_index0), levels=paste("G", 1:length(GLassoSequence)))
  

  
  if (length(LassoSequence)>=2 & length(GLassoSequence)>=2){ #### CASE1: multiple lasso and glasso
    
    lasso_index0 <- rep(LassoSequence, length(GLassoSequence))
    Glasso_index0 <- rep(1:length(GLassoSequence), each=length(LassoSequence))
    Glasso_index0 <- factor(paste("G", round(Glasso_index0, digits = 3)), levels=paste("G", 1:length(GLassoSequence)))
    
    lowestPress <- min(vec_PRESS)
    
    if(length(which(vec_PRESS == lowestPress))>1){
      #this happens when min(vec_PRESS) contains multiple values
      lowestplus1SE <- lowestPress + min(vec_se[which(vec_PRESS == lowestPress)]) 
    } else{
      lowestplus1SE <- lowestPress + vec_se[which(vec_PRESS == lowestPress)] #plot 1SE rule region 
    }
    lasso1 <- array()
    lasso2 <- array()
      for(i in 1:length(GLassoSequence)){
       
        pressindex <- which(abs(PRESS[, i]-lowestplus1SE)==min(abs(PRESS[, i]-lowestplus1SE)))  #comment: it seems that we have to have an index number here, instead of inserting the index directly in PRESS[index, i]
        pressindex <- pressindex[length(pressindex)]  #this is because in case of large Glasso values, preindex is a vector, we choose the one with the most sparse results
        lasso1temp <- LassoSequence[pressindex]
        
        if(PRESS[pressindex, i] - lowestplus1SE > 0 ){
          if(LassoSequence[pressindex] == LassoSequence[1]){
            lasso2temp <- lasso1temp  #otherwise lasso2 is out of the boundary
          } else{
            lasso2temp <- LassoSequence[pressindex - 1]
            if(PRESS[pressindex - 1, i]  - lowestplus1SE > 0){
              lasso2temp <- lasso1temp
            }
          }
          #the following condition concerns a rare case 
          
          lasso1[i] <- lasso2temp
          lasso2[i] <- lasso1temp
        }  else if (PRESS[pressindex, i] - lowestplus1SE < 0 ){
              if(LassoSequence[pressindex] == LassoSequence[length(LassoSequence)]){
                lasso2temp <- lasso1temp #otherwise lasso2 is out of the boundary
              } else{
                lasso2temp <- LassoSequence[pressindex + 1]
                #the following condition concerns a rare case 
                if(PRESS[pressindex + 1, i]  - lowestplus1SE < 0) {
                  lasso2temp <- lasso1temp
                }
              }
            
            lasso1[i] <- lasso1temp
            lasso2[i] <- lasso2temp
          
        } else { #this is when a PRESS value lies exactly on the 1SE dotted line 
      
          lasso1[i] <- lasso1temp
          lasso2[i] <- lasso1temp
        }
      }
 
    lambdaregion <- cbind(lasso1, lasso2)
    l1matrix <- t(cbind(lasso1, matrix(NA, length(lasso1), length(LassoSequence)-1)))
    l2matrix <- t(cbind(lasso2, matrix(NA, length(lasso2), length(LassoSequence)-1)))
  
    
  } else if(length(LassoSequence)>=2 & length(GLassoSequence)==1){ #### CASE2: multiple lasso, one glasso
    
    df <- data.frame(LassoI = LassoSequence, Press = vec_PRESS, Upper = upper, Lower = lower)
    
    lowestPress <- min(vec_PRESS)
    lowestplus1SE <- lowestPress + vec_se[which(vec_PRESS == lowestPress)] #plot 1SE rule region 
    lasso1 <- df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))]
    if(vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE > 0 ){
      if(df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$LassoI[1]){
        lasso2 <- lasso1 #otherwise lasso2 is out of the boundary
      } else{
        lasso2 <- df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]
      }
      #the following condition concerns a rare case 
      if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]  - lowestplus1SE > 0) &
          (df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$LassoI[1])){
        lasso2 <- lasso1
        cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", lasso1)
        lambdaregion <- glasso1
      } else{
        lambdaregion <- c(lasso2, lasso1)
      }
    } else if (vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE < 0 ){
        if(df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$LassoI[length(LassoSequence)]){
          lasso2 <- lasso1 #otherwise lasso2 is out of the boundary
        } else{
          lasso2 <- df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]
        }
      
      #the following condition concerns a rare case 
      if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]  - lowestplus1SE < 0) & 
          (df$LassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$LassoI[length(LassoSequence)])){
        lasso2 <- lasso1
        cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", lasso1)
        lambdaregion <- lasso1
      } else{
        lambdaregion <- c(lasso1, lasso2)
      }
      
    } else {#this is when a PRESS value lies exactly on the 1SE dotted line
      lasso2 <- lasso1
      lambdaregion <- lasso1
    }
    
    
  } else if(length(LassoSequence)==1 & length(GLassoSequence)>= 2){####CASE 3: Multiple glasso, one lasso
     
      df <- data.frame(GLassoI = GLassoSequence, Press = vec_PRESS, Upper = upper, Lower = lower)
    
      lowestPress <- min(vec_PRESS)
      lowestplus1SE <- lowestPress + vec_se[which(vec_PRESS == lowestPress)] #plot 1SE rule region 
      glasso1 <- df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))]
      if(vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE > 0 ){
        if(df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$GLassoI[1]){
          glasso2 <- glasso1  #otherwise Glasso2 is out of the boundary
        } else{
          glasso2 <- df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]
        }
        #the following condition concerns a rare case 
        if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) - 1]  - lowestplus1SE > 0) &
            (df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$GLassoI[1])){
          glasso2 <- glasso1
          cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", glasso1)
          lambdaregion <- glasso1
        } else{
          lambdaregion <- c(glasso2, glasso1)
        }
        
      } else if (vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] - lowestplus1SE < 0 ){
        if(df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] == df$GLassoI[length(GLassoSequence)]){
          glasso2 <- glasso1 #otherwise glasso2 is out of the boundary
        } else{
          glasso2 <- df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]
        }
        
        #the following condition concerns a rare case 
        if((vec_PRESS[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE))) + 1]  - lowestplus1SE < 0) & 
            (df$GLassoI[which(abs(vec_PRESS-lowestplus1SE)==min(abs(vec_PRESS-lowestplus1SE)))] != df$GLassoI[length(GLassoSequence)])){
          glasso2 <- glasso1
          cat("\nWarning! The region for proper tuning parameter values is not available. A possible value is ", glasso1)
          lambdaregion <- glasso1
        } else{
          lambdaregion <- c(glasso1, glasso2)
        }
        
      } else { #this is when a PRESS value lies exactly on the 1SE dotted line 
        glasso2 <- glasso1
        lambdaregion <- glasso1
      }
    
  }

  indexTuning <- which(abs(PRESS - lowestplus1SE)==min(abs(PRESS - lowestplus1SE)), arr.ind = T)
  bestTunning <- matrix(NA, dim(indexTuning)[1], 2)
  bestTunning[, 1] <- LassoSequence[indexTuning[1]]
  bestTunning[, 2] <- GLassoSequence[indexTuning[2]]
  
  colnames(bestTunning) <- c("Lasso", "Group Lasso")
  
  colnames(lambdaregion) <- c("lower bound", "upper bound")
  
  return_crossvali <- list()
  return_crossvali$PRESS <- PRESS
  return_crossvali$SE_MSE <- se_MSE
  return_crossvali$Press1SE <- lowestplus1SE
  return_crossvali$VarSelected <- varselected
  return_crossvali$plot <- p
  return_crossvali$Lasso_values <- LassoSequence
  return_crossvali$Glasso_values <- GLassoSequence
  return_crossvali$Lambdaregion <- lambdaregion
  return_crossvali$RecommendedLambda <- bestTunning
  return_crossvali$plotlog <- plotlog
  attr(return_crossvali, "class") <- "sparseSCA"
  return(return_crossvali)

}

