#'Cross-validation plot
#'
#'A plot of mean square errors + 1 standard error against 
#'Lasso tuning parameters. The plot is plotted against a log 
#'scale of lambda if \code{LassoSequence} is not defined by users. 
#'
#'@param x A object for plot.
#'@param ...  Argument to be passed to or from other methods. 
#'
#'@examples
#'\dontrun{
#'## S3 method for class 'CVstructuredSCA'
#'plot(x)
#'}
#'
#'@export
plot.CVstructuredSCA <- function(x, ...){
  
  LassoSequence <- x$LassoSequence
  PRESS <- x$MSPE
  sd_MSE <- x$Standard_Error
  plotlog <- x$plotlog
  
  upper <- PRESS + sd_MSE
  lower <- PRESS - sd_MSE
  
  LassoI <- NA
  Press <- NA
  Lower <- NA
  Upper <- NA
  
  if(plotlog == 1){
    df <- data.frame( LassoI = log(LassoSequence), Press = PRESS, Upper = upper, Lower = lower)
    lowestPress <- min(PRESS)
    lowestplus1SE <- lowestPress + min(sd_MSE[which(PRESS == lowestPress)]) #plot 1SE rule region the idea is to fine the region of lasso 
    #where according to the 1SE rule the lasso should be in that region. 
    #min() is used in case there are multiple PRESS==lowestPress
    pressindex <- which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))
    lasso1 <- df$LassoI[utils::tail(pressindex)] #in case multiple Lasso values are available, choose the most sparse one. 
    
    if(PRESS[utils::tail(pressindex)] - lowestplus1SE > 0 ){
      if(df$LassoI[utils::tail(pressindex)] == df$LassoI[1]){
        #this is the first value of the vector
        lregion <- c(lasso1, lasso1)
      } else{
        lasso2 <- df$LassoI[utils::tail(pressindex) - 1]
        lregion <- c(lasso2, lasso1)
      }
      
    } else if (PRESS[utils::tail(pressindex)] - lowestplus1SE < 0 ){
      if(df$LassoI[utils::tail(pressindex)] == df$LassoI[length(df$LassoI)]){
        #this is the last value of the vector
        lregion <- c(lasso1, lasso1)
      } else{
        lasso2 <- df$LassoI[utils::tail(pressindex) + 1]
        lregion <- c(lasso1, lasso2)
      }
      
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
    p <- p + ggplot2::labs(x = "Lasso (log scale)", y="Mean Squared Prediction Errors + 1SE")
  } else{
    
    df <- data.frame( LassoI = LassoSequence, Press = PRESS, Upper = upper, Lower = lower)
    
    lowestPress <- min(PRESS)
    lowestplus1SE <- lowestPress + min(sd_MSE[which(PRESS == lowestPress)]) #plot 1SE rule region the idea is to fine the region of lasso where according to the 1SE rule the lasso should be in that region. 
    
    pressindex <- which(abs(PRESS-lowestplus1SE)==min(abs(PRESS-lowestplus1SE)))
    lasso1 <- df$LassoI[utils::tail(pressindex)] #in case multiple lasso values available
    
    if(PRESS[utils::tail(pressindex)] - lowestplus1SE > 0 ){
      if(df$LassoI[utils::tail(pressindex)] == df$LassoI[1]){
        # this is the first element of the vector
        lregion <- c(lasso1, lasso1)
      } else{
        lasso2 <- df$LassoI[utils::tail(pressindex) - 1]
        lregion <- c(lasso2, lasso1)
      }
      
    } else if (PRESS[utils::tail(pressindex)] - lowestplus1SE < 0 ){
      if(df$LassoI[utils::tail(pressindex)] == length(LassoSequence)){
        lasso2 <- lasso1
      } else{
        lasso2 <- df$LassoI[utils::tail(pressindex) + 1]
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
    p <- p + ggplot2::labs(x = "Lasso", y="Mean Squared Prediction Errors + 1SE")
  }
  
  
  if(plotlog == 1){
    #in this case lregion is on log scale and has to be converted. 
    lregion <- exp(lregion)
  }
  return(p)
}