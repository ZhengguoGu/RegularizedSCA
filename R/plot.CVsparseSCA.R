#'Ploting Cross-validation results
#'
#'A plot of PRESS +/- 1 standard error against Lasso OR Group 
#'Lasso tuning parameters, with the vertical dotted black line 
#'indicating the lowest PRESS+1SE. 
#'
#'In case both the Lasso sequence and the Group Lasso sequence 
#'contain more than 2 elements, the cross-validation plot is replaced with 
#'a heatmap of mean squared prediction errors (MSPE) against Lasso and Group 
#'Lasso tuning parameters (x-axis: the Group Lasso; y-axis: the Lasso)
#'
#'@param x A object for plot.
#'@param ...  Argument to be passed to or from other methods. 
#'
#'@examples
#'\dontrun{
#'## S3 method for class 'CVsparseSCA'
#'plot(x)
#'}
#'
#'@export
plot.CVsparseSCA <- function(x, ...){
  
  LassoSequence <- x$Lasso_values
  GLassoSequence <- x$Glasso_values
  PRESS <- x$MSPE
  vec_PRESS <- c(PRESS)
  se_MSE <- x$SE_MSE
  vec_se <- c(se_MSE)
  
  upper <- vec_PRESS + vec_se
  lower <- vec_PRESS - vec_se 
  
  plotlog <- x$plotlog
  
  
  if (length(LassoSequence)>=2 & length(GLassoSequence)>=2){ #### CASE1: multiple lasso and glasso
    
    length_lasso <- length(x$Lasso_values)
    lasso <- 1:length_lasso
    length_glasso <- length(x$Glasso_values)
    glasso <- 1:length_glasso
    
    grid <- expand.grid(ROW=lasso, COL=glasso)
    grid$HIGHT <- c(x$MSPE)
    Acol.regions <- colorspace::diverge_hsv(n=50)
    colramp <- grDevices::colorRampPalette(c('yellow',  'green',  'blue'))
    
    if(plotlog == 1){
      p <- lattice::levelplot(HIGHT~COL*ROW, grid, col.regions = colramp,
                              contour=T, labels=T,
                              scales = list(x=list(at = c(1, length_glasso), labels=c(round(x$Glasso_values[1], digits = 3), round(x$Glasso_values[length_glasso], digits = 3))),
                                            y=list(at = c(1, length_lasso), labels=c(round(x$Lasso_values[1], digits = 3), round(log(x$Lasso_values[length_lasso]), digits = 3)))),
                              xlab = "Group Lasso tuning parameters", 
                              ylab = "Lasso tuning parameters (equal spaced on the log scale)",
                              main = "Mean squared prediction errors")
    }else{
      p <- lattice::levelplot(HIGHT~COL*ROW, grid, col.regions = colramp,
                              contour=T, labels=T,
                              scales = list(x=list(at = c(1, length_glasso), labels=c(round(x$Glasso_values[1], digits = 3), round(x$Glasso_values[length_glasso], digits = 3))),
                                            y=list(at = c(1, length_lasso), labels=c(round(x$Lasso_values[1], digits = 3), round(x$Lasso_values[length_lasso], digits = 3)))),
                              xlab = "Group Lasso tuning parameters", 
                              ylab = "Lasso tuning parameters",
                              main = "Mean Squared Prediction Errors")
      
    }
    
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
    
    if (plotlog == 1){
      df$LassoI = log(LassoSequence)
      p <- ggplot2::ggplot(df, ggplot2::aes_string(x='LassoI',y='Press')) +
        ggplot2::geom_errorbar(ggplot2::aes_string(ymin='Lower',ymax='Upper'), width=.1) +
        ggplot2::geom_point(ggplot2::aes_string(x='LassoI',y='Press'))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3)+
        #ggplot2::scale_x_discrete(limits=lasso_index[1:length(LassoSequence)]) +
        ggplot2::geom_vline(xintercept = log(lasso1), 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = log(lasso2), 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Lasso (on log scale)", y="Prediction Mean Squared Errors +/- 1SE")
    } else{
      p <- ggplot2::ggplot(df, ggplot2::aes_string(x='LassoI',y='Press')) +
        ggplot2::geom_errorbar(ggplot2::aes_string(ymin='Lower',ymax='Upper'), width=.1) +
        ggplot2::geom_point(ggplot2::aes_string(x='LassoI',y='Press'))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3)+
        #ggplot2::scale_x_discrete(limits=lasso_index[1:length(LassoSequence)]) +
        ggplot2::geom_vline(xintercept = lasso1, 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = lasso2, 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Lasso", y="Prediction Mean Squared Errors +/- 1SE")
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
    
    if (plotlog == 1){
      df$GLassoI = log(GLassoSequence)
      p <- ggplot2::ggplot(df, ggplot2::aes_string(x='GLassoI',y='Press')) +
        ggplot2::geom_errorbar(ggplot2::aes_string(ymin='Lower',ymax='Upper'), width=.1) +
        ggplot2::geom_point(ggplot2::aes_string(x='GLassoI',y='Press'))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3) +
        #ggplot2::scale_x_discrete(limits=Glasso_index[1:length(GLassoSequence)])
        ggplot2::geom_vline(xintercept = log(glasso1), 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = log(glasso2), 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Group Lasso (on log scale)", y="Prediction Mean Squared Errors +/- 1SE")
    } else{
      p <- ggplot2::ggplot(df, ggplot2::aes_string(x='GLassoI',y='Press')) +
        ggplot2::geom_errorbar(ggplot2::aes_string(ymin='Lower',ymax='Upper'), width=.1) +
        ggplot2::geom_point(ggplot2::aes_string(x='GLassoI',y='Press'))+
        ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3) +
        #ggplot2::scale_x_discrete(limits=Glasso_index[1:length(GLassoSequence)])
        ggplot2::geom_vline(xintercept = glasso1, 
                            linetype = "longdash", col = "red") +
        ggplot2::geom_vline(xintercept = glasso2, 
                            linetype = "longdash", col = "red") 
      p <- p + ggplot2::labs(x = "Group Lasso", y="Prediction Mean Squared Errors +/- 1SE")
    }
    
  }
  
  return(p)
}