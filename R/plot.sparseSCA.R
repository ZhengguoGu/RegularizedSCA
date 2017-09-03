plot.sparseSCA <- function(obj_sSCA){
  
  LassoSequence <- obj_sSCA$Lasso_values
  GLassoSequence <- obj_sSCA$Glasso_values
  PRESS <- obj_sSCA$PRESS
  vec_PRESS <- c(PRESS)
  se_MSE <- obj_sSCA$SE_MSE
  vec_se <- c(se_MSE)
  
  upper <- vec_PRESS + vec_se
  lower <- vec_PRESS - vec_se 
  vec_varsel <- c(obj_sSCA$VarSelected)
  
  plotlog <- obj_sSCA$plotlog
  

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
    
    l1s <- NA #It seems that a variable must be defined first, otherwise R CMD check generates a Note
    l2s <- NA
    
    if(plotlog == 1){
      
      df <- data.frame(GLassoI = Glasso_index0, LassoI = log(lasso_index0), Press = vec_PRESS, Upper = upper, Lower = lower, l1s = c(log(l1matrix)), l2s = c(log(l2matrix)))
      df2 <- data.frame(GLassoI = Glasso_index0, LassoI = log(lasso_index0),  Var = vec_varsel, l1s = c(log(l1matrix)), l2s = c(log(l2matrix)))
      xtag <- "log(Lasso)"
    } else{
      df <- data.frame(GLassoI = Glasso_index0, LassoI = lasso_index0, Press = vec_PRESS, Upper = upper, Lower = lower, l1s = c(l1matrix), l2s = c(l2matrix))
      df2 <- data.frame(GLassoI = Glasso_index0, LassoI = lasso_index0,  Var = vec_varsel, l1s = c(l1matrix), l2s = c(l2matrix))
      xtag <- "Lasso"
    }
    
    
    p1 <- ggplot2::ggplot(df, ggplot2::aes_string(x='LassoI',y='Press',group='GLassoI')) +
      ggplot2::facet_grid(.~GLassoI)+
      ggplot2::geom_errorbar(ggplot2::aes_string(ymin='Lower',ymax='Upper', group='GLassoI'), width=.1) +
      ggplot2::geom_point(ggplot2::aes_string(x='LassoI',y='Press',group='GLassoI')) +
      ggplot2::geom_hline(yintercept = upper[which(vec_PRESS == min(vec_PRESS))], linetype = 3) +
      ggplot2::geom_vline(data = subset(df, !is.na(l1s)), ggplot2::aes_string(xintercept = 'l1s'), linetype = "longdash", col = "red" ) +
      ggplot2::geom_vline(data = subset(df, !is.na(l2s)), ggplot2::aes_string(xintercept = 'l2s'), linetype = "longdash", col = "red" )      
    
    p1 <- p1 + ggplot2::labs(x = xtag, y="Prediction Mean Squared Errors +/- 1SE")
    
    p2 <- ggplot2::ggplot(df2, ggplot2::aes_string(x='LassoI',y='vec_varsel',group='GLassoI')) +
      ggplot2::facet_grid(.~GLassoI)+
      ggplot2::geom_point(ggplot2::aes_string(x='LassoI',y='vec_varsel',group='GLassoI')) +
      ggplot2::geom_vline(data = subset(df, !is.na(l1s)), ggplot2::aes_string(xintercept = 'l1s'), linetype = "longdash", col = "red" ) +
      ggplot2::geom_vline(data = subset(df, !is.na(l2s)), ggplot2::aes_string(xintercept = 'l2s'), linetype = "longdash", col = "red" )  
    p2 <- p2 + ggplot2::labs(x = xtag, y="# of non-zero component loadings selected in P matrix")
    
    p <- list()
    p[[1]] <- p1
    p[[2]] <- p2
    
    names(p) <- c("Cross-validation curve", "# of variables selected")
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